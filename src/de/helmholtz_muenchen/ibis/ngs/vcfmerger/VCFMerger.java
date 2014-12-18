package de.helmholtz_muenchen.ibis.ngs.vcfmerger;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.concurrent.ExecutionException;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;

import de.helmholtz_muenchen.ibis.utils.threads.Executor;
import de.helmholtz_muenchen.ibis.utils.threads.UnsuccessfulExecutionException;


public class VCFMerger {

	
		
		public static String mergeVCFs(String GATK,String RefGenome, String Infolder, String Outfolder, String Regex,String GenotypeMergeOption,final ExecutionContext exec, NodeLogger logger, String OUTFILETAG) throws InvalidSettingsException{
			
			LinkedList<String> Files2Merge = new LinkedList<String>();
			de.helmholtz_muenchen.ibis.utils.ngs.FileSearch.searchWithV(Infolder,Regex,Files2Merge);
			String OUTFILE = merge_vcfs(GATK,RefGenome, Files2Merge,Outfolder,GenotypeMergeOption,exec,logger,OUTFILETAG);
			return OUTFILE;
		}
		
		
		/**
		 * Executes vcf-merge 
		 * @param Files2Merge
		 * @param Outfolder
		 * @param exec
		 * @param logger
		 */
		private static String merge_vcfs(String GATK,String RefGenome, LinkedList<String> Files2Merge,String Outfolder,String GenotypeMergeOption, ExecutionContext exec, NodeLogger logger,String OUTFILETAG){
			ArrayList<String> command = new ArrayList<String>();

			String OUTFILE = Outfolder+"/AllSamples_vcfmerger"+OUTFILETAG+".vcf";
			String ERRFILE = Outfolder+"/AllSamples_vcfmerger"+OUTFILETAG+".vcf.err";
			
			command.add("java");
	    	command.add("-jar "+GATK);
	    	command.add("-T CombineVariants");
	    	command.add("-R "+RefGenome);
			command.addAll(Files2Merge);
			command.add("--genotypemergeoption "+GenotypeMergeOption);
			command.add("-o "+OUTFILE);
			
			try {
				Executor.executeCommand(new String[]{StringUtils.join(command, " ")}, exec,new String[]{}, logger,OUTFILE,ERRFILE);
				
			} catch (CanceledExecutionException | InterruptedException
					| ExecutionException | UnsuccessfulExecutionException e) {
				e.printStackTrace();
			}
			return OUTFILE;
		}
		
}
