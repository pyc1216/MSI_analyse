

rule all:
    input:
        'msisensor1_top100/sample.result',
        'msisensor1_pcr/sample.result',
        'msisensor1_select138/sample.result',
        'msisensor1_select138_default/sample.result',

        'msisensor2_top100/sample.result',
        'msisensor2_pcr/sample.result',
        'msisensor2_select138/sample.result',
        'msisensor2_select138_default/sample.result',

        'msisensor2_to_default/sample.result',
        'msisensor2_to/sample.result',


        'mantis_select138/sample.result',
        'mantis_top100/sample.result',
        'mantis_pcr/sample.result',


        


rule msisensor1_top100:
    input:
        t='bam/tumor.markdup.bam',
        n='bam/normal.markdup.bam',
        list='/data_sas/ypu/git_repository/HEMECDx/test_data/MSI_analyse/ucsc.hg19.microsatellites.top100.list'
    output:
        'msisensor1_top100/sample.result'
    threads:
        16
    log:
        'logs/msisensor1_top100.log'
    shell:
        """
        msisensor msi -l 1 -p 1 -q 1 -s 1 -b {threads} -d {input.list} -t {input.t} -n {input.n} -o {output} 1>{log} 2>&1
        """

rule msisensor1_pcr:
    input:
        t='bam/tumor.markdup.bam',
        n='bam/normal.markdup.bam',
        list='/data_sas/ypu/git_repository/HEMECDx/test_data/MSI_analyse/ucsc.hg19.microsatellites.pcr.list'
    output:
        'msisensor1_pcr/sample.result'
    threads:
        16
    log:
        'logs/msisensor1_pcr.log'
    shell:
        """
        msisensor msi -l 1 -p 1 -q 1 -s 1 -b {threads} -d {input.list} -t {input.t} -n {input.n} -o {output} 1>{log} 2>&1
        """
        
        
rule msisensor1_select138:
    input:
        t='bam/tumor.markdup.bam',
        n='bam/normal.markdup.bam',
        list='/data_sas/ypu/git_repository/HEMECDx/test_data/MSI_analyse/ucsc.hg19.microsatellites.select138.list'
    output:
        'msisensor1_select138/sample.result'
    threads:
        16
    log:
        'logs/msisensor1_select138.log'
    shell:
        """
        msisensor msi -l 1 -p 1 -q 1 -s 1 -b {threads} -d {input.list} -t {input.t} -n {input.n} -o {output} 1>{log} 2>&1
        """

rule msisensor1_select138_default:
    input:
        t='bam/tumor.markdup.bam',
        n='bam/normal.markdup.bam',
        list='/data_sas/ypu/git_repository/HEMECDx/test_data/MSI_analyse/ucsc.hg19.microsatellites.select138.list'
    output:
        'msisensor1_select138_default/sample.result'
    threads:
        16
    log:
        'logs/msisensor1_select138_default.log'
    shell:
        """
        msisensor msi -b {threads} -d {input.list} -t {input.t} -n {input.n} -o {output} 1>{log} 2>&1
        """

#####################################################################
rule msisensor2_top100:
    input:
        t='bam/tumor.markdup.bam',
        n='bam/normal.markdup.bam',
        list='/data_sas/ypu/git_repository/HEMECDx/test_data/MSI_analyse/ucsc.hg19.microsatellites.top100.list'
    output:
        'msisensor2_top100/sample.result'
    threads:
        16
    log:
        'logs/msisensor2_top100.log'
    shell:
        """
        /data_sas/ypu/git_repository/msisensor2/msisensor2 \
        msi -l 1 -p 1 -q 1 -s 1 -b {threads} -d {input.list} -t {input.t} -n {input.n} -o {output} 1>{log} 2>&1
        """

rule msisensor2_pcr:
    input:
        t='bam/tumor.markdup.bam',
        n='bam/normal.markdup.bam',
        list='/data_sas/ypu/git_repository/HEMECDx/test_data/MSI_analyse/ucsc.hg19.microsatellites.pcr.list'
    output:
        'msisensor2_pcr/sample.result'
    threads:
        16
    log:
        'logs/msisensor2_pcr.log'
    shell:
        """
        /data_sas/ypu/git_repository/msisensor2/msisensor2 \
        msi -l 1 -p 1 -q 1 -s 1 -b {threads} -d {input.list} -t {input.t} -n {input.n} -o {output} 1>{log} 2>&1
        """
        
rule msisensor2_select138:
    input:
        t='bam/tumor.markdup.bam',
        n='bam/normal.markdup.bam',
        list='/data_sas/ypu/git_repository/HEMECDx/test_data/MSI_analyse/ucsc.hg19.microsatellites.select138.list'
    output:
        'msisensor2_select138/sample.result'
    threads:
        16
    log:
        'logs/msisensor2_select138.log'
    shell:
        """
        /data_sas/ypu/git_repository/msisensor2/msisensor2 \
        msi -l 1 -p 1 -q 1 -s 1 -b {threads} -d {input.list} -t {input.t} -n {input.n} -o {output} 1>{log} 2>&1
        """
        
rule msisensor2_select138_default:
    input:
        t='bam/tumor.markdup.bam',
        n='bam/normal.markdup.bam',
        list='/data_sas/ypu/git_repository/HEMECDx/test_data/MSI_analyse/ucsc.hg19.microsatellites.select138.list'
    output:
        'msisensor2_select138_default/sample.result'
    threads:
        16
    log:
        'logs/msisensor2_select138_default.log'
    shell:
        """
        /data_sas/ypu/git_repository/msisensor2/msisensor2 \
        msi -b {threads} -d {input.list} -t {input.t} -n {input.n} -o {output} 1>{log} 2>&1
        """
        
        
        
        
#####################################################################
rule msisensor2_to_default:
    input:
        t='bam/tumor.markdup.bam',
        model='/data_sas/ypu/git_repository/msisensor2/models_hg19_GRCh37'
    output:
        'msisensor2_to_default/sample.result'
    threads:
        16
    log:
        'logs/msisensor2_to_default.log'
    shell:
        """
        /data_sas/ypu/git_repository/msisensor2/msisensor2 \
        msi -b {threads} -t {input.t} -o {output} -M {input.model} 1>{log} 2>&1
        """

rule msisensor2_to:
    input:
        t='bam/tumor.markdup.bam',
        model='/data_sas/ypu/git_repository/msisensor2/models_hg19_GRCh37',
    output:
        'msisensor2_to/sample.result'
    threads:
        16
    log:
        'logs/msisensor2_to.log'
    shell:
        """
        /data_sas/ypu/git_repository/msisensor2/msisensor2 \
        msi -l 1 -p 1 -q 1 -s 1 -b {threads} -t {input.t} -o {output} -M {input.model} 1>{log} 2>&1
        """

#####################################################################
rule mantis_select138:
    input:
        t='bam/tumor.markdup.bam',
        n='bam/normal.markdup.bam',
    output:
        'mantis_select138/sample.result'
    threads:
        16
    log:
        'logs/mantis_select138.log'
    shell:
        """
        python /data_sas/ypu/git_repository/MANTIS/mantis.py \
        --bedfile /data_sas/ypu/git_repository/HEMECDx/test_data/MSI_analyse/med1cdx_ms.bed \
        --genome /home/ypu//4-Ref/hg19/ucsc.hg19.fasta \
        --threads {threads} \
        -n {input.n} \
        -t {input.t} \
        -o {output} 1>{log} 2>&1
        """

rule mantis_top100:
    input:
        t='bam/tumor.markdup.bam',
        n='bam/normal.markdup.bam',
    output:
        'mantis_top100/sample.result'
    threads:
        16
    log:
        'logs/mantis_top100.log'
    shell:
        """
        python /data_sas/ypu/git_repository/MANTIS/mantis.py \
        --bedfile /data_sas/ypu/git_repository/HEMECDx/test_data/MSI_analyse/mantis.top100.bed \
        --genome /home/ypu//4-Ref/hg19/ucsc.hg19.fasta \
        --threads {threads} \
        -n {input.n} \
        -t {input.t} \
        -o {output} 1>{log} 2>&1
        """

rule mantis_pcr:
    input:
        t='bam/tumor.markdup.bam',
        n='bam/normal.markdup.bam',
    output:
        'mantis_pcr/sample.result'
    threads:
        16
    log:
        'logs/mantis_pcr.log'
    shell:
        """
        python /data_sas/ypu/git_repository/MANTIS/mantis.py \
        --bedfile /data_sas/ypu/git_repository/HEMECDx/test_data/MSI_analyse/mantis.pcr.bed\
        --genome /home/ypu//4-Ref/hg19/ucsc.hg19.fasta \
        --threads {threads} \
        -n {input.n} \
        -t {input.t} \
        -o {output} 1>{log} 2>&1
        """
#####################################################################

