# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.8
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.





from sys import version_info
if version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_Chemical', [dirname(__file__)])
        except ImportError:
            import _Chemical
            return _Chemical
        if fp is not None:
            try:
                _mod = imp.load_module('_Chemical', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _Chemical = swig_import_helper()
    del swig_import_helper
else:
    import _Chemical
del version_info
try:
    _swig_property = property
except NameError:
    pass  # Python < 2.2 doesn't have 'property'.


def _swig_setattr_nondynamic(self, class_type, name, value, static=1):
    if (name == "thisown"):
        return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name, None)
    if method:
        return method(self, value)
    if (not static):
        if _newclass:
            object.__setattr__(self, name, value)
        else:
            self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)


def _swig_setattr(self, class_type, name, value):
    return _swig_setattr_nondynamic(self, class_type, name, value, 0)


def _swig_getattr_nondynamic(self, class_type, name, static=1):
    if (name == "thisown"):
        return self.this.own()
    method = class_type.__swig_getmethods__.get(name, None)
    if method:
        return method(self)
    if (not static):
        return object.__getattr__(self, name)
    else:
        raise AttributeError(name)

def _swig_getattr(self, class_type, name):
    return _swig_getattr_nondynamic(self, class_type, name, 0)


def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except Exception:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object:
        pass
    _newclass = 0


class Chemical(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Chemical, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Chemical, name)
    __repr__ = _swig_repr
    __swig_setmethods__["m_mw"] = _Chemical.Chemical_m_mw_set
    __swig_getmethods__["m_mw"] = _Chemical.Chemical_m_mw_get
    if _newclass:
        m_mw = _swig_property(_Chemical.Chemical_m_mw_get, _Chemical.Chemical_m_mw_set)
    __swig_setmethods__["m_K_ow"] = _Chemical.Chemical_m_K_ow_set
    __swig_getmethods__["m_K_ow"] = _Chemical.Chemical_m_K_ow_get
    if _newclass:
        m_K_ow = _swig_property(_Chemical.Chemical_m_K_ow_get, _Chemical.Chemical_m_K_ow_set)
    __swig_setmethods__["m_pKa"] = _Chemical.Chemical_m_pKa_set
    __swig_getmethods__["m_pKa"] = _Chemical.Chemical_m_pKa_get
    if _newclass:
        m_pKa = _swig_property(_Chemical.Chemical_m_pKa_get, _Chemical.Chemical_m_pKa_set)
    __swig_setmethods__["m_frac_non_ion"] = _Chemical.Chemical_m_frac_non_ion_set
    __swig_getmethods__["m_frac_non_ion"] = _Chemical.Chemical_m_frac_non_ion_get
    if _newclass:
        m_frac_non_ion = _swig_property(_Chemical.Chemical_m_frac_non_ion_get, _Chemical.Chemical_m_frac_non_ion_set)
    __swig_setmethods__["m_frac_unbound"] = _Chemical.Chemical_m_frac_unbound_set
    __swig_getmethods__["m_frac_unbound"] = _Chemical.Chemical_m_frac_unbound_get
    if _newclass:
        m_frac_unbound = _swig_property(_Chemical.Chemical_m_frac_unbound_get, _Chemical.Chemical_m_frac_unbound_set)
    __swig_setmethods__["m_r_s"] = _Chemical.Chemical_m_r_s_set
    __swig_getmethods__["m_r_s"] = _Chemical.Chemical_m_r_s_get
    if _newclass:
        m_r_s = _swig_property(_Chemical.Chemical_m_r_s_get, _Chemical.Chemical_m_r_s_set)
    __swig_setmethods__["m_acid_base"] = _Chemical.Chemical_m_acid_base_set
    __swig_getmethods__["m_acid_base"] = _Chemical.Chemical_m_acid_base_get
    if _newclass:
        m_acid_base = _swig_property(_Chemical.Chemical_m_acid_base_get, _Chemical.Chemical_m_acid_base_set)

    def __init__(self):
        this = _Chemical.new_Chemical()
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_destroy__ = _Chemical.delete_Chemical
    __del__ = lambda self: None

    def Init(self, arg2, arg3, arg4, arg5, arg6, arg7):
        return _Chemical.Chemical_Init(self, arg2, arg3, arg4, arg5, arg6, arg7)

    def InitConfig(self, conf):
        return _Chemical.Chemical_InitConfig(self, conf)

    def calcIon(self):
        return _Chemical.Chemical_calcIon(self)

    def calcBinding(self):
        return _Chemical.Chemical_calcBinding(self)
Chemical_swigregister = _Chemical.Chemical_swigregister
Chemical_swigregister(Chemical)

# This file is compatible with both classic and new-style classes.


