package de.helmholtz_muenchen.ibis.utils;

import java.io.File;

public class BinaryHandler {

	public final static String TOOL_PATH = "/home/software/bin";
	
	@Deprecated
	public static String checkToolAvailability(String BinaryName) {
		return checkToolAvailability(BinaryName, TOOL_PATH);
	}
	
	public static String checkToolAvailability(String BinaryName, String dir){
		
		FileSearch f = new FileSearch();
		f.searchDirectory(new File(dir), BinaryName);
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
