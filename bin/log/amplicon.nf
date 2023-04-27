import Logos

logo = new Logo()
println '\n'
println logo.show()

def param_log(){
log.info """
AMPLICON PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                      ${params.workflow}
LOG IS TBD

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________
"""

}
