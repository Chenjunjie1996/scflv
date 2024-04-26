import groovy.json.JsonSlurper


class Helper {
    public static Map getProtocol(workflow, log, protocol) {
        def jsonSlurper = new JsonSlurper()
        def json = new File("${workflow.projectDir}/assets/protocols.json").text
        def protocols = jsonSlurper.parseText(json)
        if(protocols.containsKey(protocol)) {
            return protocols[protocol]
        } else {
            log.warn("Protocol '${protocol}' not recognized by the pipeline. Passing on the protocol to the aligner unmodified.")
            return ["protocol": protocol]
        }
    }
}
