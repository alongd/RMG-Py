#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2026 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

"""Discover reaction-direction, reverse-rate, and reaction I/O sites.

Uses Python tokens rather than importing code, so Cython .pyx sources receive
exactly the same coverage as Python. Comments and ordinary strings are excluded;
executable Jinja template expressions are scanned separately.
Site identities include the enclosing function, primitive, and occurrence number;
line numbers are evidence only and do not invalidate verdicts after unrelated edits.
"""
import argparse
import ast
from collections import Counter
import io
import json
from pathlib import Path
import re
import tokenize

DIRECTION = {'is_isomorphic', 'is_same_reaction', 'matches_species', 'has_template',
             'check_identical', 'same_species_lists', 'are_identical_species_references',
             'allows_reverse_match', 'find_degenerate_reactions', 'check_for_duplicates',
             'mark_valid_duplicates', 'mark_duplicate_reaction', 'mark_duplicate_reactions',
             'chemkin_duplicate_flags', 'add_path_reaction', 'combine_models',
             'compare_model_reactions', 'log_enlarge_summary', 'add_reverse_attribute',
             'convert_duplicates_to_multi', 'reaction_matches', 'get_reaction_matches'}
REVERSE = {'check_reaction_rate_serialization', 'check_resolved_species_reversibility',
           'check_reverse_from_equilibrium_supported',
           'calculate_microcanonical_rate_coefficient', 'calculate_rate_coefficients',
           'get_rate_coefficients_CSE_Advanced', 'get_rate_coefficients_SLS',
           'generate_kinetics', 'fit_interpolation_model', 'fit_interpolation_models',
           'get_equilibrium_constants', '_reverse_rate_at_reference_potential',
           'get_equilibrium_constant', 'generate_reverse_rate_coefficient',
           'reverse_arrhenius_rate', 'reverse_surface_arrhenius_rate',
           'reverse_sticking_coeff_rate', 'reverse_surface_charge_transfer_rate',
           'reverse_arrhenius_charge_transfer_rate'}
SERIALIZE = {'to_cantera', 'to_chemkin', 'write_kinetics_entry', 'read_kinetics_entry',
             'read_reaction_comments', 'load_chemkin_file', 'render_chemkin_file',
             'save_chemkin_file', 'save_chemkin_surface_file', 'reaction_to_dicts',
             'reaction_to_dict_list', 'get_reaction_equation', '_build_equation_string',
             'generate_cantera_data', 'get_mech_dict', 'obj_to_dict', 'as_dict', 'make_object',
             'recursive_make_object', 'expand_to_dict', 'reaction_state',
             'apply_reaction_state', 'copy_reaction'}
BUILD = {'Reaction', 'LibraryReaction', 'TemplateReaction', 'DepositoryReaction',
         'PDepReaction', 'make_new_reaction'}
DEFINITION = re.compile(r'^( *)(?:(?:async )?def|cpdef|cdef)\s+(?:(?:[\w.*]+)\s+)?(\w+)\s*\(')
CLASS = re.compile(r'^( *)(?:cdef )?class\s+(\w+)')


IO_PREFIXES = ('read_', 'write_', 'load_', 'save_', 'render_', 'export',
               'serialize', 'deserialize', 'parse_')
IO_NAMES = {'save', 'load', 'write', 'read', '__repr__', '__reduce__',
            'as_dict', 'make_object', 'get_mech_dict', 'generate_cantera_data'}
GENERIC_IO = {'as_dict', 'make_object', 'recursive_make_object', 'expand_to_dict',
              'reaction_state', 'apply_reaction_state', 'copy_reaction'}
REACTION_IDENTIFIERS = {'reaction', 'reactions', 'reaction_list', 'net_reactions',
                        'path_reactions', 'reaction_model', 'reaction_systems',
                        'Reaction', 'LibraryReaction', 'TemplateReaction',
                        'reactants', 'products', 'kinetics'}
PRIMITIVES = DIRECTION | REVERSE | SERIALIZE | BUILD
REVERSE_NAME = re.compile(r'(?:generate_reverse|reverse_).*?(?:rate|kinetics).*')

ARKANE_SERIALIZATION_CHOKEPOINT = 'check_reaction_rate_serialization'
ARKANE_EMITTER_CALLS = {
    'save_kinetics_lib', 'save_yaml', 'to_chemkin', 'write_kinetics_entry',
}
ARKANE_EMITTER_FUNCTIONS = {
    'save', 'save_input_file', 'save_kinetics_lib', 'write_chemkin', 'write_output',
}
ARKANE_REACTION_NAMES = {
    'kinetics', 'net_reactions', 'path_reactions', 'reaction', 'reaction_str',
    'reactions', 'rxn', 'rxn_list',
}
ARKANE_RATE_NAMES = {'kinetics', 'rate', 'rates', 'reaction', 'rxn'}
ARKANE_WRITE_CALLS = {'dump', 'save', 'write', 'writelines'}
ARKANE_SYNTAX_MARKERS = ('kinetics(', 'pdepreaction(', 'reaction(')


def is_reaction_io(owner, body):
    """Include named I/O or declaration-bearing functions and generic serializers."""
    name = owner.rsplit('.', 1)[-1]
    names = {token.string for token in body}
    declaration = 'reversible' in names and bool(names & {'reactants', 'products', 'kinetics'})
    io = name.startswith(IO_PREFIXES) or name in SERIALIZE | IO_NAMES or declaration
    return io and (name in GENERIC_IO or bool(names & REACTION_IDENTIFIERS))


def discover(root):
    """Return every primitive call, including calls nested in Cython functions."""
    result = []
    for package in ('rmgpy', 'arkane'):
        for path in sorted((root / package).rglob('*')):
            if path.suffix not in ('.py', '.pyx'):
                continue
            source = path.read_text()
            tokens = list(tokenize.generate_tokens(io.StringIO(source).readline))
            # Enclosing definitions are taken from tokenized code lines, excluding
            # docstrings, comments, and continuation lines inside expressions.
            definitions = {}
            stack = []
            for token in tokens:
                if token.type != tokenize.NAME:
                    continue
                line = token.line.rstrip('\n')
                prefix = line[:token.start[1]]
                if prefix.strip():
                    continue
                indent = len(prefix)
                while stack and indent <= stack[-1][0]:
                    stack.pop()
                match = DEFINITION.match(line) or CLASS.match(line)
                if match:
                    stack.append((indent, match.group(2)))
                definitions[token.start[0]] = '.'.join(name for _, name in stack) or '<module>'
            scope = '<module>'
            counts = Counter()
            significant = [t for t in tokens if t.type not in (
                tokenize.ENCODING, tokenize.NL, tokenize.NEWLINE, tokenize.INDENT,
                tokenize.DEDENT, tokenize.COMMENT, tokenize.ENDMARKER)]
            for i, token in enumerate(significant):
                scope = definitions.get(token.start[0], scope)
                name = token.string
                if (token.type != tokenize.NAME or i + 1 == len(significant)
                        or significant[i + 1].string != '('):
                    continue
                if name not in PRIMITIVES and not REVERSE_NAME.fullmatch(name):
                    continue
                if i and significant[i-1].string in ('def', 'cpdef', 'cdef', 'class'):
                    continue
                # Typed Cython definitions also are declarations, not calls.
                if (DEFINITION.match(token.line) and token.start[1] <= token.line.index('(')):
                    continue
                kind = ('direction' if name in DIRECTION else 'reverse' if name in REVERSE
                        else 'builder' if name in BUILD else 'reader-exporter' if name in SERIALIZE else 'reverse')
                key = (scope, name)
                counts[key] += 1
                identity = f'{path.relative_to(root)}::{scope}::{name}#{counts[key]}'
                result.append({'site': identity, 'kind': kind, 'line': token.start[0],
                               'source': token.line.strip()})
            # Jinja expressions in HTML templates are executable source too.
            # Inspect expressions, while continuing to exclude prose/docstrings.
            for token in tokens:
                if token.type != tokenize.STRING or not any(marker in token.string for marker in ('{{', '{%')):
                    continue
                owner = definitions.get(token.start[0], '<module>')
                for expression in re.findall(r'\{[{%](.*?)[}%]\}', token.string, re.S):
                    for name in sorted(DIRECTION | REVERSE | SERIALIZE | BUILD):
                        for match in re.finditer(r'\b' + name + r'\s*\(', expression):
                            counts[(owner, name)] += 1
                            result.append({'site': f'{path.relative_to(root)}::{owner}::{name}#{counts[(owner, name)]}',
                                           'kind': 'template-expression',
                                           'line': token.start[0],
                                           'source': expression.strip()})
            # Reader/exporter entry points with reaction-bearing identifiers are
            # inventoried even when they delegate serialization and never construct
            # a Reaction directly (e.g. pickle/repr and library save_entry).
            functions = {}
            current = '<module>'
            for token in tokens:
                current = definitions.get(token.start[0], current)
                if token.type == tokenize.NAME:
                    functions.setdefault(current, []).append(token)
            for owner, body in functions.items():
                if not is_reaction_io(owner, body):
                    continue
                token = body[0]
                result.append({'site': f'{path.relative_to(root)}::{owner}::reaction-io',
                               'kind': 'reader-exporter', 'line': token.start[0],
                               'source': token.line.strip()})
    return sorted(result, key=lambda row: row['site'])


def _function_scopes(tree):
    """Yield each Python function with its class/function-qualified owner."""
    result = []

    def visit(body, parents=()):
        for node in body:
            if isinstance(node, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
                owner = parents + (node.name,)
                if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    result.append(('.'.join(owner), node))
                visit(node.body, owner)

    visit(tree.body)
    return result


def discover_arkane_emitters(root):
    """Discover Arkane functions that directly emit or delegate reaction/rate text."""
    result = []
    for path in sorted((root / 'arkane').rglob('*.py')):
        tree = ast.parse(path.read_text(), filename=str(path))
        for owner, function in _function_scopes(tree):
            names = {node.id for node in ast.walk(function) if isinstance(node, ast.Name)}
            attributes = {node.attr for node in ast.walk(function) if isinstance(node, ast.Attribute)}
            calls = [node for node in ast.walk(function) if isinstance(node, ast.Call)]
            call_names = []
            for call in calls:
                if isinstance(call.func, ast.Name):
                    call_names.append(call.func.id)
                elif isinstance(call.func, ast.Attribute):
                    call_names.append(call.func.attr)
            strings = [node.value for node in ast.walk(function)
                       if isinstance(node, ast.Constant) and isinstance(node.value, str)]
            signals = set()
            for marker in ARKANE_SYNTAX_MARKERS:
                if any(marker in value for value in strings):
                    signals.add('syntax:' + marker)
            for call_name in ARKANE_EMITTER_CALLS.intersection(call_names):
                signals.add('call:' + call_name)
            basename = owner.rsplit('.', 1)[-1]
            identifiers = names | attributes
            if (basename in ARKANE_EMITTER_FUNCTIONS
                    and ARKANE_WRITE_CALLS.intersection(call_names)
                    and ARKANE_REACTION_NAMES.intersection(identifiers)
                    and ARKANE_RATE_NAMES.intersection(identifiers)):
                signals.add('reaction-rate-writer')
            if basename == 'save_kinetics_lib':
                signals.add('kinetics-library-writer')
            if owner == 'Arkane.get_libraries':
                signals.add('explorer-object-emitter')
            if not signals:
                continue
            chokepoint_calls = [call for call in calls
                                if ((isinstance(call.func, ast.Name)
                                     and call.func.id == ARKANE_SERIALIZATION_CHOKEPOINT)
                                    or (isinstance(call.func, ast.Attribute)
                                        and call.func.attr == ARKANE_SERIALIZATION_CHOKEPOINT))]
            emission_calls = [
                call for call, call_name in zip(calls, call_names)
                if call_name in {'open', 'save', 'write', 'write_kinetics_entry'}
            ]
            chokepoint_checks = []
            for call in chokepoint_calls:
                arguments = [ast.unparse(argument) for argument in call.args]
                arguments.extend(
                    f'{keyword.arg}={ast.unparse(keyword.value)}' for keyword in call.keywords)
                chokepoint_checks.append(arguments)
            result.append({
                'site': f'{path.relative_to(root)}::{owner}',
                'line': function.lineno,
                'signals': sorted(signals),
                'chokepoint_calls': len(chokepoint_calls),
                'chokepoint_arities': sorted(len(call.args) + len(call.keywords)
                                             for call in chokepoint_calls),
                'chokepoint_checks': sorted(chokepoint_checks),
                'chokepoint_lines': sorted(call.lineno for call in chokepoint_calls),
                'first_emission_line': min((call.lineno for call in emission_calls), default=None),
            })
    return sorted(result, key=lambda row: row['site'])


def validate_arkane_emitters(emitters, inventory):
    """Keep the discovered Arkane writer inventory and chokepoint routing complete."""
    discovered = {row['site'] for row in emitters}
    recorded = set(inventory)
    errors = []
    if discovered - recorded:
        errors.append('Unclassified Arkane emitters: ' + ', '.join(sorted(discovered - recorded)))
    if recorded - discovered:
        errors.append('Disappeared Arkane emitters: ' + ', '.join(sorted(recorded - discovered)))
    for site, decision in inventory.items():
        if decision.get('verdict') not in ('checked', 'excluded'):
            errors.append('Invalid Arkane emitter verdict: ' + site)
        if not decision.get('reason') or not decision.get('tests'):
            errors.append('Missing Arkane emitter reason or test: ' + site)
    for row in emitters:
        decision = inventory.get(row['site'], {})
        if decision.get('verdict') == 'checked':
            if row['chokepoint_calls'] < 1:
                errors.append('Emitter bypasses serialization chokepoint: ' + row['site'])
            elif any(arity != 3 for arity in row['chokepoint_arities']):
                errors.append('Emitter omits exact rate or direction at chokepoint: ' + row['site'])
            elif row['chokepoint_checks'] != sorted(decision.get('checks', [])):
                errors.append('Emitter checks the wrong rate or direction at chokepoint: ' + row['site'])
            elif (row['first_emission_line'] is not None
                  and max(row['chokepoint_lines']) >= row['first_emission_line']):
                errors.append('Emitter reaches serialization before its chokepoint: ' + row['site'])
    if errors:
        raise ValueError('\n'.join(errors))


def validate(sites, verdicts):
    """Reject unclassified new sites, stale verdicts, and incomplete decisions."""
    discovered = {row['site'] for row in sites}
    recorded = set(verdicts)
    errors = []
    if discovered - recorded:
        errors.append('Unclassified sites: ' + ', '.join(sorted(discovered - recorded)))
    if recorded - discovered:
        errors.append('Disappeared sites: ' + ', '.join(sorted(recorded - discovered)))
    for site, verdict in verdicts.items():
        if verdict.get('verdict') not in ('guarded', 'refused', 'not-reachable-for-resolved'):
            errors.append('Invalid verdict: ' + site)
        if not verdict.get('reason') or not verdict.get('tests'):
            errors.append('Missing reason or test: ' + site)
    if errors:
        raise ValueError('\n'.join(errors))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--verdicts', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    sites = discover(args.root)
    if args.verdicts:
        verdicts = json.loads(args.verdicts.read_text())
        validate(sites, verdicts)
        sites = [dict(row, **verdicts[row['site']]) for row in sites]
    payload = json.dumps(sites, indent=2) + '\n'
    if args.output:
        args.output.write_text(payload)
    else:
        print(payload, end='')


if __name__ == '__main__':
    main()
