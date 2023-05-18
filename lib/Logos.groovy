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

class Logo {
    def c = new Colors()

    def logColors() {
        // Mod Log colors if not ANSI
        final mono = !env.NXF_ANSI_LOG
        //final mono = true
        if ( mono ) {
            for (k in c.keySet()) {
                c[k] = ''
            }
        }
    }

    private def frameLogo(String logo, String color=c.blu) {
        // return framed lines of passed logo text
        def logos = logo.split('[\n\r]')
        int maxLen = ( logos.collect{ it.length() } ).max() + 4
        def frameit = { it * maxLen }
        // def bar = "${c.dim}.${c.non}"
        def bar = ""
        def logoLines = ''
        for ( line in logos ) {
            (line =~ /\S/) && \
            (logoLines += "${bar}${color} ${line}   ${c.non}${bar}\n")
        }
        logos = "${c.dim}.${frameit('.')}.${c.non}\n" \
              + logoLines \
              + "${bar}${c.dim}${frameit('.')}${bar}"
        logos = logoLines
        return logos.stripIndent()
    }

    public def show(String logo=this.logo, String clr=c.non) {
        // return framed lines of chosen org logo in chosen color
        frameLogo(logo, clr)
    }

/*
ASCII Name Art Options: http://patorjk.com/software/taag/
*/

static def logo = this.logo_jaxngsops_cyber

static def logo_jaxgm_ansi_regular = $/
     ██  █████  ██   ██        ██████ ███████       ███    ██  ██████  ███████        ██████  ██████  ███████ 
     ██ ██   ██  ██ ██        ██      ██            ████   ██ ██       ██            ██    ██ ██   ██ ██      
     ██ ███████   ███   █████ ██      ███████ █████ ██ ██  ██ ██   ███ ███████ █████ ██    ██ ██████  ███████ 
██   ██ ██   ██  ██ ██        ██           ██       ██  ██ ██ ██    ██      ██       ██    ██ ██           ██ 
 █████  ██   ██ ██   ██        ██████ ███████       ██   ████  ██████  ███████        ██████  ██      ███████ 
/$

static def logo_jaxcsngsops_big = $/
     _   _   __  __      ____ ____        _   _  ____ ____         ___  ____  ____  
    | | / \  \ \/ /     / ___/ ___|      | \ | |/ ___/ ___|       / _ \|  _ \/ ___| 
 _  | |/ _ \  \  /_____| |   \___ \ _____|  \| | |  _\___ \ _____| | | | |_) \___ \ 
| |_| / ___ \ /  |_____| |___ ___) |_____| |\  | |_| |___) |_____| |_| |  __/ ___) |
 \___/_/   \_/_/\_\     \____|____/      |_| \_|\____|____/       \___/|_|   |____/ 
/$

static def logo_jaxngsops_mini = $/
                 _  __         __  __    _   _   __ 
   |  /\  \/ __ /  (_ __ |\ | /__ (_ __ / \ |_) (_  
 \_| /--\ /\    \_ __)   | \| \_| __)   \_/ |   __) 
/$


static def logo_jaxngsops_cyber = $/
_____ _______ _     _     _______ _______     __   _  ______ _______      _____   _____  _______
  |   |_____|  \___/  ___ |       |______ ___ | \  | |  ____ |______ ___ |     | |_____| |______
__|   |     | _/   \_     |_____  ______|     |  \_| |_____| ______|     |_____| |       ______|
/$                                                                                               



}
