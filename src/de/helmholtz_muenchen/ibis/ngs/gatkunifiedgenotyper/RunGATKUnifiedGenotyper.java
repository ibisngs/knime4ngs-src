package de.helmholtz_muenchen.ibis.ngs.gatkunifiedgenotyper;


import java.io.File;

import org.knime.core.node.ExecutionContext;

import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.threads.Executor;


public class RunGATKUnifiedGenotyper extends GATKUnifiedGenotyperNodeModel{
	
protected void CallVariants(ExecutionContext exec, String gatk, String[] bam, String ref, String out, String intf, String dbsnpf, String snpIndel, int threads, double [] param, String baq, boolean filter, String proxyOptions, int GATK_MEMORY_USAGE, String opt_flags) throws Exception {
		
		String lockFile = out + SuccessfulRunChecker.LOCK_ENDING;
	
		//for each thread 2G
		String cmd="java -jar -Xmx"+GATK_MEMORY_USAGE+"G " + proxyOptions + gatk;
		cmd+=" -T UnifiedGenotyper";
		cmd+=" -nt "+threads;
		
		for(int i=0;i<bam.length;i++){
			cmd+=" -I "+bam[i];
		}
		
		cmd+=" -R "+ref;
		cmd+=" -o "+out;
		cmd+=" -glm "+snpIndel;
		
		if(!dbsnpf.equals("")){
			cmd+=" -D "+dbsnpf;
		}
		
		if(!intf.equals("")){
			cmd+=" -L "+intf;
		}
		
		cmd+=" -stand_call_conf "+param[0];
		cmd+=" -stand_emit_conf "+param[1];
		cmd+=" -pcr_error "+param[2];
		cmd+=" -contamination "+param[3];
		cmd+=" -hets "+param[4];
		cmd+=" -deletions "+param[5];
		cmd+=" -mbq "+(int) param[6];
		cmd+=" -indelHeterozygosity "+param[7];
		cmd+=" -minIndelCnt "+ (int) param[8]; 
		cmd+=" -minIndelFrac "+param[9];
		cmd+=" -indelGOP "+ (int) param[10];
		cmd+=" -indelGCP "+(int) param[11];
		
		cmd+=" -baq "+baq;
		
		if(filter){
			cmd+=" --filter_mismatching_base_and_quals ";
		}
		
		if(!opt_flags.equals("")) {
			cmd += " "+opt_flags;
		}
		
		// run command
		GATKUnifiedGenotyperNodeModel.logger.info("Running GATK UnifiedGenotyper...");
		GATKUnifiedGenotyperNodeModel.logger.info("Log files can be found in "+out+".out.log and "+out+".err.log");
		
//		Executor.executeCommand(new String[]{cmd}, exec, null, GATKUnifiedGenotyperNodeModel.logger, out+".out.log", out+".err.log", null);
		super.executeCommand(new String[]{cmd}, exec, new File(lockFile),out+".out.log", out+".err.log");
	}

}
