package de.helmholtz_muenchen.ibis.utils;

import java.io.File;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;

public class BinaryHandler {

	
	public static String checkToolAvailability(String BinaryName){
		
		String TOOL_DIR = IBISKNIMENodesPlugin.getDefault().getToolDirPreference();
		
		FileSearch f = new FileSearch();
		f.searchDirectory(new File(TOOL_DIR), BinaryName);
		f.getResult();

		String ToolPath = null;
		int count = f.getResult().size();
		if(count ==0){
		    System.out.println("\nNo result found!");
		}else{
		    System.out.println("\nFound " + count + " result!\n");
		    for (String matched : f.getResult()){
		    	ToolPath = matched;
		    	System.out.println("Found : " + matched);
		    }
		}
		
		return ToolPath;
	}
	
	
	
	
	
	
}
