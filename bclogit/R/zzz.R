if (getRversion() >= "2.15.1") {
    utils::globalVariables(c("_bclogit_process_matched_pairs_cpp"))
}

.onAttach <- function(libname, pkgname) {
    packageStartupMessage(
        paste("Welcome to bclogit v",
            utils::packageVersion("bclogit"),
            sep = ""
        )
    )
}
