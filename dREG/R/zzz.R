.First.lib <-
    function(libname, pkgname, where)
    library.dynam("dREG", pkgname, libname)

.Last.lib <-
    function(libpath)
    dyn.unload(file.path(libpath,
                         "libs",
                         paste("dREG",
                               .Platform$"dynlib.ext",
                               sep = "")))
