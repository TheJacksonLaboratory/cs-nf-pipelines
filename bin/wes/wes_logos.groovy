#!/usr/bin/env groovy
//
// Copyright 2019 Benjamin Leopold
//          _    _
//   _ _ _ | |_ |_| ___  _ _ _  ___  ___  ___  ___  ___ 
//  | ' ' || , \| |/ , \| ' ' |/ -_|/  / / , \|  _|/ -_|
//  |_|_|_||___/|_|\___/|_|_|_|\___.\__\ \___/|_|  \___.
//


//NB! relies on Nextflow env variable NXF_ANSI_LOG
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

    private def frameLogo(String logo, String color=c.grn) {
        // return framed lines of passed logo text
        def logos = logo.split('[\n\r]')
        int maxLen = ( logos.collect{ it.length() } ).max() + 4
        def frameit = { it * maxLen }
        def bar = "${c.dim}|${c.non}"
        def logoLines = ''
        for ( line in logos ) {
            (line =~ /\S/) && \
            (logoLines += "${bar}${color} ${line}   ${c.non}${bar}\n")
        }
        logos = "${c.dim},${frameit('_')},${c.non}\n" \
              + logoLines \
              + "${bar}${c.dim}${frameit('_')}${bar}"
        return logos.stripIndent()
    }

    public def show(String logo=this.logo, String clr=c.grn) {
        // return framed lines of chosen org logo in chosen color
        frameLogo(logo, clr)
    }

/*
ASCII Name Art Options:

( Many produced using: http://patorjk.com/software/taag/ )
*/

static def logo = this.logo_jaxgm_block

static def logo_jaxgm_block = $/
                                                 _       _ 
       _|    _|_|    _|      _|          _|_|_|  _|      _|
       _|  _|    _|    _|  _|          _|        _|_|  _|_|
       _|  _|_|_|_|      _|    _|_|_|  _|  _|_|  _|  _|  _|
 _|    _|  _|    _|    _|  _|          _|    _|  _|      _|
   _|_|    _|    _|  _|      _|          _|_|_|  _|      _|
/$

static def logo_jaxgm_banner3 = $/
                                                     
      ##    ###    ##     ##      ######   ##     ## 
      ##   ## ##    ##   ##      ##    ##  ###   ### 
      ##  ##   ##    ## ##       ##        #### #### 
      ## ##     ##    ###   #### ##   #### ## ### ## 
##    ## #########   ## ##       ##    ##  ##     ## 
##    ## ##     ##  ##   ##      ##    ##  ##     ## 
 ######  ##     ## ##     ##      ######   ##     ## 
/$

static def logo_jaxgm_arrows = $/
                                                                  
     >=>       >>       >=>      >=>       >===>    >=>       >=> 
     >=>      >>=>       >=>   >=>       >>    >=>  >> >=>   >>=> 
     >=>     >> >=>       >=> >=>       >=>         >=> >=> > >=> 
     >=>    >=>  >=>        >=>    >==> >=>         >=>  >=>  >=> 
     >=>   >=====>>=>     >=> >=>       >=>   >===> >=>   >>  >=> 
>>   >=>  >=>      >=>   >=>   >=>       >=>    >>  >=>       >=> 
 >===>   >=>        >=> >=>      >=>      >====>    >=>       >=> 
/$

static def logo_jaxgm_diamond = $/
                                                                  
     /\\      /\       /\\      /\\         /\\\\   /\\       /\\ 
     /\\     /\ \\      /\\   /\\         /\    /\\ /\ /\\   /\\\ 
     /\\    /\  /\\      /\\ /\\         /\\        /\\ /\\ / /\\ 
     /\\   /\\   /\\       /\\    /\\\\\ /\\        /\\  /\\  /\\ 
     /\\  /\\\\\\ /\\    /\\ /\\         /\\   /\\\ /\\   /\  /\\ 
/\   /\\ /\\       /\\  /\\   /\\         /\\   /\  /\\       /\\ 
 /\\\\  /\\         /\\/\\      /\\        /\\\\\   /\\       /\\ 
/$


static def logo_mbcore1 = $/
         _    _                                       
  _ _ _ | |_ |_| ___  _ _ _  ___  ___  ___  ___  ___  
 | ' ' || . \| |/ . \| ' ' |/ ._|/  / / . \|  _|/ ._| 
 |_|_|_||___/|_|\___/|_|_|_|\___.\__\ \___/|_|  \___. 
/$

static def logo_mbcore2 = $/
       __     __         ___  __   __   __   ___ 
 |\/| |__) | /  \  |\/| |__  /  ` /  \ |__) |__  
 |  | |__) | \__/  |  | |___ \__, \__/ |  \ |___ 
/$

static def logo_mbcore3 = $/
                                         
  ._ _  |_  o  _  ._ _   _   _  _  ._ _  
  | | | |_) | (_) | | | (/_ (_ (_) | (/_ 
/$

static def logo_mbcore4 = $/
          __                                                                    
         /\ \      __                                                           
  ___ ___\ \ \____/\_\    ___     ___ ___      __    ___    ___   _ __    __    
/' __` __`\ \ '__`\/\ \  / __`\ /' __` __`\  /'__`\ /'___\ / __`\/\`'__\/'__`\  
/\ \/\ \/\ \ \ \L\ \ \ \/\ \L\ \/\ \/\ \/\ \/\  __//\ \__//\ \L\ \ \ \//\  __/  
\ \_\ \_\ \_\ \_,__/\ \_\ \____/\ \_\ \_\ \_\ \____\ \____\ \____/\ \_\\ \____\ 
 \/_/\/_/\/_/\/___/  \/_/\/___/  \/_/\/_/\/_/\/____/\/____/\/___/  \/_/ \/____/ 
/$

static def logo_mbcore5 = $/
           _     _                                         
 _ __ ___ | |__ (_) ___  _ __ ___   ___  ___ ___  _ __ ___ 
| '_ ` _ \| '_ \| |/ _ \| '_ ` _ \ / _ \/ __/ _ \| '__/ _ \
| | | | | | |_) | | (_) | | | | | |  __/ (_| (_) | | |  __/
|_| |_| |_|_.__/|_|\___/|_| |_| |_|\___|\___\___/|_|  \___|
/$

}


// vim: set ft=groovy ts=4 sw=0 tw=100 et fdm=syntax:
