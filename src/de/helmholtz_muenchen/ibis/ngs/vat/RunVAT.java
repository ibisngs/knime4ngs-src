package de.helmholtz_muenchen.ibis.ngs.vat;


import org.knime.core.node.ExecutionContext;

import de.helmholtz_muenchen.ibis.utils.threads.Executor;


public class RunVAT {
	
	public static void IndelMapper(ExecutionContext exec, String indelmapper, String intervals, String transcripts, String vcfin, String vcfout) throws Exception{
		
		String cmd=indelmapper;
		cmd+=" "+intervals;
		cmd+=" "+transcripts;
		
		VATNodeModel.logger.info("Running VAT indelMapper...");
		VATNodeModel.logger.info("Log file can be found in "+vcfout+".err.log");
		
		Executor.executeCommand(new String []{cmd}, exec, null, VATNodeModel.logger, vcfout, vcfout+".err.log", vcfin);
		
	}
	
	public static void SNPMapper(ExecutionContext exec, String snpmapper, String intervals, String transcripts, String vcfin, String vcfout)throws Exception{
		
		String cmd=snpmapper;
		cmd+=" "+intervals;
		cmd+=" "+transcripts;
		
		VATNodeModel.logger.info("Running VAT snpMapper...");
		VATNodeModel.logger.info("Log file can be found in "+vcfout+".err.log");
		
		Executor.executeCommand(new String []{cmd}, exec, null, VATNodeModel.logger, vcfout, vcfout+".err.log", vcfin);

	}

}