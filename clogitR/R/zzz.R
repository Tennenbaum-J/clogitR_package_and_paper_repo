.onAttach = function(libname, pkgname){
	packageStartupMessage(
			paste("Welcome to clogitR v", 
					utils::packageVersion("clogitR"), 
					sep = "")
	)	
}