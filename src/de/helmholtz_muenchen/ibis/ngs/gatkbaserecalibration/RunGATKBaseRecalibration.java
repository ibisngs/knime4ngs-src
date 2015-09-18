package de.helmholtz_muenchen.ibis.ngs.gatkbaserecalibration;


import java.io.File;
import java.util.Map;

import org.knime.core.node.ExecutionContext;

import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;

public class RunGATKBaseRecalibration extends GATKBaseRecalibrationNodeModel  {
	
	// first call of base recalibrator
	protected void BaseRecalibrator(ExecutionContext exec, String gatk, String inputbam, String inputref, String outtable, String pphase1, String pmills, String pdbsnp, String pint, boolean [] cov, int tailqual, double go, int maxcycle, int [] indelmis, int cputhreads, String proxyOptions, int GATK_MEMORY_USAGE, String opt_flags) throws Exception {
		
		String lockFile = outtable + SuccessfulRunChecker.LOCK_ENDING;
		
		//create command string
		String cmd="java -jar -Xmx"+GATK_MEMORY_USAGE+"G "+ proxyOptions + gatk;
		cmd+=" -T BaseRecalibrator";
		cmd+=" -R "+inputref;
		cmd+=" -I "+inputbam;
		cmd+=" -o "+outtable;
		
		if(!pphase1.equals("")){
			cmd+=" -knownSites "+pphase1;
		}
		
		if(!pmills.equals("")){
			cmd+=" -knownSites "+pmills;
		}
		
		if(!pdbsnp.equals("")){
			cmd+=" -knownSites "+pdbsnp;
		}
		
		if(!pint.equals("")){
			cmd+=" -L "+pint;
		}
		
		cmd+=" -noStandard";
		if(cov[0]){
			cmd+=" -cov ContextCovariate";
		}
		if(cov[1]){
			cmd+=" -cov CycleCovariate";
		}
		if(cov[2]){
			cmd+=" -cov QualityScoreCovariate";
		}
		if(cov[3]){
			cmd+=" -cov ReadGroupCovariate";
		}
		if(cov[4]){
			cmd+=" -cov RepeatLengthCovariate";
		}
		if(cov[5]){
			cmd+=" -cov RepeatUnitCovariate";
		}
		
		cmd+=" -lqt "+tailqual;
		cmd+=" -bqsrBAQGOP "+go;
		cmd+=" -maxCycle "+maxcycle;
		
		cmd+=" -ddq "+indelmis[0];
		cmd+=" -idq "+indelmis[1];
		cmd+=" -mdq "+indelmis[2];
		cmd+=" -ics "+indelmis[3];
		cmd+=" -mcs "+indelmis[4];
		
		cmd+=" -nct "+cputhreads;
		
		if(!opt_flags.equals("")) {
			cmd+=" "+opt_flags;
		}
		
		// run command
		GATKBaseRecalibrationNodeModel.logger.info("Running GATK BaseRecalibrator...");
		GATKBaseRecalibrationNodeModel.logger.info("Log files can be found in "+outtable+".out.log and "+outtable+".err.log");
		
		GATKBaseRecalibrationNodeModel.logger.info(cmd);
		
//		Executor.executeCommand(new String [] {cmd} ,exec, null, GATKBaseRecalibrationNodeModel.logger, outtable+".out.log", outtable+".err.log", null);
		super.executeCommand(new String[]{cmd}, exec, new File(lockFile),outtable+".out.log", outtable+".err.log");
	}
	
	
	// second call of base recalibrator
	protected void BaseRecalibrator(ExecutionContext exec, String gatk, String inputbam, String inputref, String inputtable, String outtable, String pphase1, String pmills, String pdbsnp, String pint, boolean [] cov, int tailqual, double go, int maxcycle, int [] indelmis, int cputhreads, String proxyOptions, int GATK_MEMORY_USAGE, String opt_flags) throws Exception {
		
		String lockFile = outtable + SuccessfulRunChecker.LOCK_ENDING;
		
		//create command string
		String cmd="java -jar -Xmx"+GATK_MEMORY_USAGE+"G "+ proxyOptions + gatk;
		cmd+=" -T BaseRecalibrator";
		cmd+=" -R "+inputref;
		cmd+=" -I "+inputbam;
		cmd+=" -BQSR "+inputtable;
		cmd+=" -o "+outtable;
		
		if(!pphase1.equals("")){
			cmd+=" -knownSites "+pphase1;
		}
		
		if(!pmills.equals("")){
			cmd+=" -knownSites "+pmills;
		}
		
		if(!pdbsnp.equals("")){
			cmd+=" -knownSites "+pdbsnp;
		}
		
		if(!pint.equals("")){
			cmd+=" -L "+pint;
		}
		
		cmd+=" -noStandard";
		if(cov[0]){
			cmd+=" -cov ContextCovariate";
		}
		if(cov[1]){
			cmd+=" -cov CycleCovariate";
		}
		if(cov[2]){
			cmd+=" -cov QualityScoreCovariate";
		}
		if(cov[3]){
			cmd+=" -cov ReadGroupCovariate";
		}
		if(cov[4]){
			cmd+=" -cov RepeatLengthCovariate";
		}
		if(cov[5]){
			cmd+=" -cov RepeatUnitCovariate";
		}
		
		cmd+=" -lqt "+tailqual;
		cmd+=" -bqsrBAQGOP "+go;
		cmd+=" -maxCycle "+maxcycle;
		
		cmd+=" -ddq "+indelmis[0];
		cmd+=" -idq "+indelmis[1];
		cmd+=" -mdq "+indelmis[2];
		cmd+=" -ics "+indelmis[3];
		cmd+=" -mcs "+indelmis[4];
		
		cmd+=" -nct "+cputhreads;
		
		cmd+= " " + opt_flags;
		// run command
		GATKBaseRecalibrationNodeModel.logger.info("Running GATK BaseRecalibrator...");
		GATKBaseRecalibrationNodeModel.logger.info("Log files can be found in "+outtable+".out.log and "+outtable+".err.log");
		
//		Executor.executeCommand(new String[]{cmd},exec, null, GATKBaseRecalibrationNodeModel.logger, outtable+".out.log", outtable+".err.log", null);
		super.executeCommand(new String[]{cmd}, exec, new File(lockFile),outtable+".out.log", outtable+".err.log");
	}
	
	protected void PrintReads(ExecutionContext exec, String gatk, String inputbam, String inputref, String inputtable, String outbam, boolean outsimple, int cputhreads, String proxyOptions, int GATK_MEMORY_USAGE) throws Exception {
		
		String lockFile = outbam + SuccessfulRunChecker.LOCK_ENDING;
		
		//create command string
		String cmd="java -jar -Xmx"+GATK_MEMORY_USAGE+"G "+ proxyOptions + gatk;
		cmd+=" -T PrintReads";
		cmd+=" -R "+inputref;
		cmd+=" -I "+inputbam;
		cmd+=" -BQSR "+inputtable;
		cmd+=" -o "+outbam;
		
		cmd+=" -nct "+cputhreads;
		
		if(outsimple){
			cmd+=" -s ";
		}
		
		//run command
		GATKBaseRecalibrationNodeModel.logger.info("Running GATK PrintReads...");
		GATKBaseRecalibrationNodeModel.logger.info("Log files can be found in "+outbam+".out.log and "+outbam+".err.log");
		
//		Executor.executeCommand(new String[]{cmd}, exec, null, GATKBaseRecalibrationNodeModel.logger, outbam+".out.log", outbam+".err.log", null);
		super.executeCommand(new String[]{cmd}, exec, new File(lockFile),outbam+".out.log", outbam+".err.log");
	}
	
	protected void AnalyzeCovariates(ExecutionContext exec, String gatk, String inputref, String beforetable, String aftertable, String pplots, String pint, String proxyOptions, int GATK_MEMORY_USAGE,String recalintermediate) throws Exception {
		
		String lockFile = pplots + SuccessfulRunChecker.LOCK_ENDING;
		
		String cmd="java -jar -Xmx"+GATK_MEMORY_USAGE+"G "+ proxyOptions + gatk;
		cmd+=" -T AnalyzeCovariates";
		cmd+=" -R "+inputref;
		cmd+=" -before "+beforetable;
		cmd+=" -after "+aftertable;
		cmd+=" -plots "+pplots;
		cmd+=" -csv "+recalintermediate;
		
		if(!pint.equals("")){
			cmd+=" -L "+pint;
		}
		
		// PATH environment is needed for calling Rscript
		Map<String, String> map =System.getenv();
		System.out.println(map.get("PATH"));
		String [] env = new String []{"PATH="+map.get("PATH")};
		
		// run command
		GATKBaseRecalibrationNodeModel.logger.info("Running GATK AnalyzeCovariates...");
		GATKBaseRecalibrationNodeModel.logger.info("Log files can be found in "+pplots+".out.log and "+pplots+".err.log");
		
//		Executor.executeCommand(new String[]{cmd}, exec, env, GATKBaseRecalibrationNodeModel.logger, pplots+".out.log", pplots+".err.log", null);
//		super.executeCommand(new String[]{cmd}, exec, new File(lockFile),pplots+".out.log", pplots+".err.log");
		super.executeCommand(new String[]{cmd}, exec, env, new File(lockFile), pplots+".out.log", pplots+".err.log", null, null, null);
	}
}
