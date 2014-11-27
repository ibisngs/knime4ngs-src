package de.helmholtz_muenchen.ibis.ngs.vcfmerger;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.concurrent.ExecutionException;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;

import de.helmholtz_muenchen.ibis.utils.threads.Executor;
import de.helmholtz_muenchen.ibis.utils.threads.UnsuccessfulExecutionException;

public class VCFMerger {

	
		
		public static void mergeVCFs(String Infolder, String Outfolder, String Regex,final ExecutionContext exec){
			
			LinkedList<String> Files2Merge = new LinkedList<String>();
			search(Infolder,Regex,Files2Merge,exec);
//			System.out.println(Arrays.toString(Files2Merge.toArray()));
		}
	
	
		/**
		 * Finds all Files in a given Directory that end with regex
		 * @param dataDirectory
		 * @param regex
		 * @param Files2Merge
		 */
		private static void search(String dataDirectory, String regex, LinkedList<String> Files2Merge,final ExecutionContext exec){
			File root = new File( dataDirectory );
		    File[] list = root.listFiles();
		        for ( File f : list ) {
		            if ( f.isDirectory() ) {
		            	search( f.getAbsolutePath() , regex, Files2Merge,exec);
		            }else {
		            	String filePath = f.getAbsolutePath();
		            	if(filePath.endsWith(regex)){			//Add if file ends with specified regex
		            		Files2Merge.add(BGZipandTabix(filePath, exec));//, new String[]{filePath,"0",null,null,"-1","0"}); 
		            	}else{
		            	}
		            }
		        } 
		}
		
		
		private static String BGZipandTabix(String Filepath,final ExecutionContext exec){
			try {
				ArrayList<String> command = new ArrayList<String>();
				command.add("/bin/sh");
				command.add("-c");
				command.add("bgzip -c "+Filepath);//+" > "+Filepath+".gz; tabix -p vcf "+Filepath+".gz");
				Executor.executeCommand(new String[]{StringUtils.join(command, " ")}, exec, null);
//				System.out.println(StringUtils.join(command, " "));
				
//				String com1="bgzip "+Filepath;
//				Process p1 = Runtime.getRuntime().exec(new String[]{"/bin/sh","-c",com1});
//				p1.waitFor();
//				System.out.println(com1);
//				Process p2 = Runtime.getRuntime().exec(new String[]{"/bin/sh","-c","tabix -p vcf "+Filepath+".gz"});
//				p2.waitFor();
//				System.out.println("tabix -p vcf "+Filepath+".gz");
				return Filepath+".gz";
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (CanceledExecutionException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (ExecutionException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (UnsuccessfulExecutionException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			return null;
		}
		
		
		
}
