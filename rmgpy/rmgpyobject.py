#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2018 Prof. William H. Green (whgreen@mit.edu),           #
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

import logging

import numpy as np

################################################################################


class RMGPYObject(object):
    """
    This class provides a general `as_dict` method to help with yml file construction
    for other objects in rmgpy to inherit
    """
    def __init__(self):
        pass

    def as_dict(self):
        """
        A helper function for YAML dumping non-cytonized objects
        """
        output_dict = dict()
        output_dict['class'] = self.__class__.__name__
        for key in self.__dict__:
            val = getattr(self, key)
            if val is not None and not callable(val) and not key.startswith('_') and val != '':
                output_dict[key] = val
        for key, val in output_dict.iteritems():
            if isinstance(val, list) and isinstance(val[0], RMGPYObject):
                output_dict[key] = [v.as_dict() for v in val]
            elif not isinstance(val, (int, float, str, list, dict)):
                if isinstance(val, np.ndarray):
                    output_dict[key] = val.tolist()
                elif val is not None:
                    # this is an object, call as_dict() again
                    output_dict[key] = val.as_dict()
        return output_dict

    def make_object(self, data, class_dict):
        """
        A helper function for YAML parsing of non-cytonized objects
        """
        for key, val in data.iteritems():
            if isinstance(val, dict) and 'class' in val:
                # Call make_object to make another object within the parent object
                class_name = val['class']
                del val['class']
                try:
                    class_to_make = class_dict[class_name]
                except KeyError:
                    raise KeyError("Class {0} must be provided in the 'class_dict' parameter "
                                   "to make the object.".format(class_name))
                obj = class_to_make()
                obj.make_object(val, class_dict)
                logging.debug("made object {0}".format(class_name))
                data[key] = obj
                # if class_name == 'TransportData':
                    # print "key: ", key
                    # print "obj: ", obj
            elif isinstance(val, list) and isinstance(val[0], dict) and 'class' in val[0]:
                # Call make_object to make a list of objects within the parent object (as in Conformer.Modes)
                data[key] = list()
                for entry in val:
                    class_name = entry['class']
                    del entry['class']
                    try:
                        class_to_make = class_dict[class_name]
                    except KeyError:
                        raise KeyError("Class {0} must be provided in the 'class_dict' parameter "
                                       "to make the object.".format(class_name))
                    obj = class_to_make()
                    obj.make_object(entry, class_dict)
                    logging.debug("made object {0}".format(class_name))
                    data[key].append(obj)
            elif isinstance(val, str):
                try:
                    float(val)
                except ValueError:
                    pass
        self.__init__(**data)
