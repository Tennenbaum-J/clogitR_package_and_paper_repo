.onAttach <- function(libname, pkgname) {
    packageStartupMessage(
        paste("Welcome to bclogit v",
            utils::packageVersion("bclogit"),
            sep = ""
        )
    )
}
