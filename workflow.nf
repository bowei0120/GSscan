#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

 def get_sample_information(sample_information_file) {
     Channel
         .fromPath(sample_information_file)
         .splitCsv(header:true, quote:'\"')
         .map{ row ->
             [row.sample_id, row.bam]
         }
 }


include { GENERATE_READ_COUNT } from "./main"
include { SINGLE_SAMPLE_GC } from "./main"
include { MERGE_ALL_SAMPLE } from "./main"


workflow{

   dat=get_sample_information(params.input)
   dat.view()

   GENERATE_READ_COUNT(dat)
   generate_rc_result=GENERATE_READ_COUNT.out.result


   SINGLE_SAMPLE_GC(generate_rc_result)

   sigle_gc_collect=SINGLE_SAMPLE_GC.out.result.collect()

   MERGE_ALL_SAMPLE(sigle_gc_collect)





}