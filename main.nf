#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/circleseq
========================================================================================
    Github : https://github.com/nf-core/circleseq
    Website: https://nf-co.re/circleseq
    Slack  : https://nfcore.slack.com/channels/circleseq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/

params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { CIRCLESEQ } from './workflows/circleseq'

//
// WORKFLOW: Run main nf-core/circleseq analysis pipeline
//
workflow NFCORE_CIRCLESEQ {
    CIRCLESEQ ()
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_CIRCLESEQ ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
