package de.helmholtz_muenchen.ibis.ngs.vcfmerger;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.concurrent.ExecutionException;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.NodeLogger;

import de.helmholtz_muenchen.ibis.utils.threads.Executor;
import de.helmholtz_muenchen.ibis.utils.threads.UnsuccessfulExecutionException;

public class VCFMerger {

	
		
		public static String mergeVCFs(String GATK,String RefGenome, String Infolder, String Outfolder, String Regex,String GenotypeMergeOption,final ExecutionContext exec, NodeLogger logger){
			
			LinkedList<String> Files2Merge = new LinkedList<String>();
			search(Infolder,Regex,Files2Merge,exec,logger);
			String OUTFILE = merge_vcfs(GATK,RefGenome, Files2Merge,Outfolder,GenotypeMergeOption,exec,logger);
			
			return OUTFILE;
		}
	
	
		
		/**
		 * Finds all Files in a given Directory that end with regex
		 * @param dataDirectory
		 * @param regex
		 * @param Files2Merge
		 * @param logger 
		 */
		private static void search(String dataDirectory, String regex, LinkedList<String> Files2Merge,final ExecutionContext exec, NodeLogger logger){
			File root = new File( dataDirectory );
		    File[] list = root.listFiles();
		        for ( File f : list ) {
		            if ( f.isDirectory() ) {
		            	search( f.getAbsolutePath() , regex, Files2Merge,exec,logger);
		            }else {
		            	String filePath = f.getAbsolutePath();
		            	if(filePath.endsWith(regex)){			//Add if file ends with specified regex
		            		Files2Merge.add("-V "+filePath);//, new String[]{filePath,"0",null,null,"-1","0"}); 
		            	}
		            }
		        } 
		}
		
	
		
		/**
		 * Executes vcf-merge 
		 * @param Files2Merge
		 * @param Outfolder
		 * @param exec
		 * @param logger
		 */
		private static String merge_vcfs(String GATK,String RefGenome, LinkedList<String> Files2Merge,String Outfolder,String GenotypeMergeOption, ExecutionContext exec, NodeLogger logger){
			ArrayList<String> command = new ArrayList<String>();

			String OUTFILE = Outfolder+"/AllSamples.vcfmerger.vcf";
			String ERRFILE = Outfolder+"/AllSamples.vcfmerger.vcf.err";
			
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
