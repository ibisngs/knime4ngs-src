package de.helmholtz_muenchen.ibis.ngs.gatkrealignment;


import org.knime.core.node.ExecutionContext;

import de.helmholtz_muenchen.ibis.utils.threads.Executor;


public class RunGATKRealignment {
	
	
	//calls gatk targetcreator
	protected static void targetcreator(ExecutionContext exec, String outputint, String inputbam, String ref, String gatk, String phase1, String mills, String interval, int threads, int maxint, int minreads, double mismatch, int window, String proxyOptions, int GATK_MEMORY_USAGE) throws Exception{
		
		//create command string
		
		// each data thread requires 2 GB of memory
		int memory = GATK_MEMORY_USAGE*threads;

		String cmd="java -jar -Xmx"+memory+"G " + proxyOptions + gatk;
		cmd+=" -T RealignerTargetCreator";
		cmd+=" -nt "+threads;
		cmd+=" -R "+ref;
		cmd+=" -I "+inputbam;
		cmd+=" -o "+outputint;
		
		if(!phase1.equals("")){
			cmd+=" -known "+phase1;
		}
		
		if(!mills.equals("")){
			cmd+=" -known "+mills;
		}
		
		if(!interval.equals("")){
			cmd+=" -L "+interval;
		}
		
		cmd+=" -maxInterval "+maxint;
		cmd+=" -minReads "+minreads;
		cmd+=" -mismatch "+mismatch;
		cmd+=" -window "+window;
		
		
		GATKRealignmentNodeModel.logger.info("Running GATK TargetCreator...");
		GATKRealignmentNodeModel.logger.info("Log files can be found in "+outputint+".out.log and "+outputint+".err.log");
		
		Executor.executeCommand(new String[]{cmd}, exec, null, GATKRealignmentNodeModel.logger, outputint+".out.log", outputint+".err.log", null);
	}
	
	//calls gatk indel realigner
	protected static void realign (ExecutionContext exec, String outputint, String outputbam, String inputbam, String ref, String gatk, String phase1, String mills, String interval, String consmode, double lod, double entropy, int maxcons, int maxisize, int maxposmove, int maxreadscons, int maxreadsaln, boolean notag, String proxyOptions, int GATK_MEMORY_USAGE, String opt_flags) throws Exception{
    	
		//create command string

		String cmd="java -jar -Xmx"+GATK_MEMORY_USAGE+"G " + proxyOptions + gatk;
		cmd+=" -T IndelRealigner";
		cmd+=" -R "+ref;
		cmd+=" -I "+inputbam;
		cmd+=" -o "+outputbam;
		cmd+=" -targetIntervals "+outputint;
		
		if(!phase1.equals("")){
			cmd+=" -known "+phase1;
		}
		
		if(!mills.equals("")){
			cmd+=" -known "+mills;
		}
		
		if(!interval.equals("")){
			cmd+=" -L "+interval;
		}
		
		cmd+=" -model "+consmode;
		cmd+=" -LOD "+lod;
		cmd+=" -entropy "+entropy;
		cmd+=" -maxConsensuses "+maxcons;
		cmd+=" -maxIsize "+maxisize;
		cmd+=" -maxPosMove "+maxposmove;
		cmd+=" -greedy "+maxreadscons;
		cmd+=" -maxReads "+maxreadsaln;
		
		if(notag){
			cmd+=" -noTags ";
		}
		
		if(!opt_flags.equals("")) {
			cmd+= " " + opt_flags;
		}
		
		
		GATKRealignmentNodeModel.logger.info("Running GATK IndelRealigner...");
		GATKRealignmentNodeModel.logger.info("Log files can be found in "+outputbam+".out.log and "+outputbam+".err.log");
		
		Executor.executeCommand(new String[] {cmd}, exec, null, GATKRealignmentNodeModel.logger, outputbam+".out.log", outputbam+".err.log", null);
	
		
	}

}
