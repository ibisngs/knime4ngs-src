package de.helmholtz_muenchen.ibis.ngs.gatkbaserecalibration;


import java.util.Map;

import org.knime.core.node.ExecutionContext;

import de.helmholtz_muenchen.ibis.utils.threads.Executor;

public class RunGATKBaseRecalibration {
	
	// first call of base recalibrator
	protected static void BaseRecalibrator(ExecutionContext exec, String gatk, String inputbam, String inputref, String outtable, String pphase1, String pmills, String pdbsnp, String pint, boolean [] cov, int tailqual, double go, int maxcycle, int [] indelmis, int cputhreads) throws Exception {
		
		//create command string
		String cmd="java -jar -Xmx4G "+gatk;
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
		
		// run command
		GATKBaseRecalibrationNodeModel.logger.info("Running GATK BaseRecalibrator...");
		GATKBaseRecalibrationNodeModel.logger.info("Log files can be found in "+outtable+".out.log and "+outtable+".err.log");
		
		GATKBaseRecalibrationNodeModel.logger.info(cmd);
		
		Executor.executeCommand(new String [] {cmd} ,exec, null, GATKBaseRecalibrationNodeModel.logger, outtable+".out.log", outtable+".err.log", null);
	}
	
	
	// second call of base recalibrator
	protected static void BaseRecalibrator(ExecutionContext exec, String gatk, String inputbam, String inputref, String inputtable, String outtable, String pphase1, String pmills, String pdbsnp, String pint, boolean [] cov, int tailqual, double go, int maxcycle, int [] indelmis, int cputhreads) throws Exception {
		
		//create command string
		String cmd="java -jar -Xmx4G "+gatk;
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
		
		
		// run command
		GATKBaseRecalibrationNodeModel.logger.info("Running GATK BaseRecalibrator...");
		GATKBaseRecalibrationNodeModel.logger.info("Log files can be found in "+outtable+".out.log and "+outtable+".err.log");
		
		Executor.executeCommand(new String[]{cmd},exec, null, GATKBaseRecalibrationNodeModel.logger, outtable+".out.log", outtable+".err.log", null);
	}
	
	protected static void PrintReads(ExecutionContext exec, String gatk, String inputbam, String inputref, String inputtable, String outbam, boolean outsimple, int cputhreads) throws Exception {
		
		//create command string
		String cmd="java -jar -Xmx4G "+gatk;
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
		
		Executor.executeCommand(new String[]{cmd}, exec, null, GATKBaseRecalibrationNodeModel.logger, outbam+".out.log", outbam+".err.log", null);
	}
	
	protected static void AnalyzeCovariates(ExecutionContext exec, String gatk, String inputref, String beforetable, String aftertable, String pplots, String pint) throws Exception {
		
		String cmd="java -jar -Xmx4G "+gatk;
		cmd+=" -T AnalyzeCovariates";
		cmd+=" -R "+inputref;
		cmd+=" -before "+beforetable;
		cmd+=" -after "+aftertable;
		cmd+=" -plots "+pplots;
		
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
		
		Executor.executeCommand(new String[]{cmd}, exec, env, GATKBaseRecalibrationNodeModel.logger, pplots+".out.log", pplots+".err.log", null);
	}
}
