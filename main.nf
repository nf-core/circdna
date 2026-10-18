#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/circdna
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/circdna
    Website: https://nf-co.re/circdna
    Slack  : https://nfcore.slack.com/channels/circdna
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

include { validateParameters; paramsHelp; paramsSummaryMap } from 'plugin/nf-schema'
include { CIRCDNA } from './workflows/circdna'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main nf-core/circdna analysis pipeline
//
workflow NFCORE_CIRCDNA {
    //   This is an example of how to use getGenomeAttribute() to fetch parameters
    //   from igenomes.config using `--genome`
    params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')

    // Print help message if needed
    if (params.help) {
        def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
        def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
        def String command = "nextflow run ${workflow.manifest.name} --input samplesheet.csv --genome GRCh38 --outdir results"
        log.info logo + paramsHelp(command: command, hideWarning: true) + citation + NfcoreTemplate.dashedLine(params.monochrome_logs)
        System.exit(0)
    }

    // Validate input parameters
    if (params.validate_params) {
        validateParameters()
    }

    WorkflowMain.initialise(workflow, params, log, args)

    CIRCDNA ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_CIRCDNA ()

    //
    // COMPLETION EMAIL AND SUMMARY
    // 关键：把隐式对象先捕获为闭包外部的局部变量。
    // 因为 onComplete 闭包在 entry workflow 退出后才执行，此时直接访问 `workflow` / `params`
    // 隐式对象会变为 null（导致 `Cannot get property 'email' on null object`）。
    // 通过 Groovy 闭包词法作用域捕获 wf / ps 引用即可绕过该问题。
    //
    def wf = workflow
    def ps = params
    def lg = log
    def pd = projectDir

    workflow.onComplete {
        def summary_params = [:]
        try {
            summary_params = paramsSummaryMap(wf)
        } catch (Exception e) {
            lg.warn "Could not parse parameters summary: ${e.message}"
        }

        try {
            if (ps.email || ps.email_on_fail) {
                NfcoreTemplate.email(wf, ps, summary_params, pd, lg, [])
            }
            NfcoreTemplate.dump_parameters(wf, ps)
            NfcoreTemplate.summary(wf, ps, lg)
            if (ps.hook_url) {
                NfcoreTemplate.IM_notification(wf, ps, summary_params, pd, lg)
            }
        } catch (Exception e) {
            lg.warn "Could not execute completion handler: ${e.message}"
        }
    }

    workflow.onError {
        try {
            if (wf.errorReport && wf.errorReport.contains("Process requirement exceeds available memory")) {
                println("🛑 Default resources exceed availability 🛑 ")
                println("💡 See here on how to configure pipeline: https://nf-co.re/docs/usage/configuration#tuning-workflow-resources 💡")
            }
        } catch (Exception e) {
            lg.warn "Could not check error report: ${e.message}"
        }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
