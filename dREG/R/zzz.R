#.First.lib and .Last.lib are obsolete and will not be used in R >= 3.0.0

#.First.lib <-
#    function(libname, pkgname, where)
#    library.dynam("dREG", pkgname, libname)

#.Last.lib <-
#    function(libpath)
#    dyn.unload(file.path(libpath,
#                         "libs",
#                         paste("dREG",
#                               .Platform$"dynlib.ext",
#                               sep = "")))

.onAttach<- function(libname, pkgName)
{
	err_cmds <- c();
	path <- Sys.which('bedmap');
	if( path == "")
	{
		#packageStartupMessage("* Can't find 'bedmap' comand which is required in some functions, please use Sys.which('bedmap') to locate it or install it.(http://bedops.readthedocs.org/en/latest/index.html)")
		err_cmds <- c(err_cmds, "bedmap");
	}

	if(length(err_cmds)>0)
	{
		packageStartupMessage("WARNING! The following dependencies were not found in your current environment.");
		packageStartupMessage("------");
		packageStartupMessage(paste(err_cmds, collapse=","));
		packageStartupMessage("------");
		packageStartupMessage("These dependencies are required by dREG to improve the processing speed.");
		packageStartupMessage("To troubleshoot these please do the following");
		packageStartupMessage("1. Make sure these commands are installed.");
		packageStartupMessage("2. check the $PATH variable to call Sys.getenv('PATH')");
		packageStartupMessage("   - or - check the dependency location by the command 'which', e.g. Sys.which('ls')");
		packageStartupMessage("3. If the command path can not be found in the default setting, plase call Sys.setenv() to set your command path.");
		packageStartupMessage("   e.g. Sys.setenv( PATH = paste(Sys.getenv('PATH'), '/your/command/path', sep=':') )" );
	}
}