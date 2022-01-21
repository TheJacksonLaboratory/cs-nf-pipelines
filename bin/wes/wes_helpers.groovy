#!/usr/bin/env groovy
//
// Copyright 2019 Benjamin Leopold
//          _    _
//   _ _ _ | |_ |_| ___  _ _ _  ___  ___  ___  ___  ___
//  | ' ' || , \| |/ , \| ' ' |/ -_|/  / / , \|  _|/ -_|
//  |_|_|_||___/|_|\___/|_|_|_|\___.\__\ \___/|_|  \___.
//


// ANSI color codes
class Colors {
    final non = "\033[0m"
    final dim = "\033[2m"
    final blk = "\033[0;30m"
    final grn = "\033[0;32m"
    final ylw = "\033[0;33m"
    final blu = "\033[0;34m"
    final pur = "\033[0;35m"
    final cyn = "\033[0;36m"
    final wht = "\033[0;37m"
}


class Summary {
    static def show(Map summary) {
        def c = new Colors()
        int len = ( summary.collect{ k,v -> k.length() } ).max()
        def frameit = { it * 79 }
        return summary.collect { k,v -> "${c.cyn}${k.padRight(len)}${c.non} : $v" }.join("\n") \
               + '\n' + frameit('_')
    }
}


class CheckParams {
    // Usage: check_numeric(Map params, List vars_require_number)
    static def is_numeric(Map params, List vars) {
        def err_numb = 'parameter must be a number. Provided value:'
        for (var in vars) {
            try {
                params[var] as Double
            } catch (e) {
                println "'${var}' ${err_numb} '${params[var]}'"
            }
        }
    }

}

