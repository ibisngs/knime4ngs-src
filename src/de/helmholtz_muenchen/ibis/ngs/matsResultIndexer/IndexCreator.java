package de.helmholtz_muenchen.ibis.ngs.matsResultIndexer;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Creates a GFF3 file which can be indexed by MISO
 * @author Michael Kluge
 *
 */
public class IndexCreator {
	
	private static final String BASE_NAME 		= "fromGTF";
	private static final String NOVE_NAME 		= ".novelEvents";
	
	private static final String FILENAME_SE 	= ".SE.txt";
	private static final String FILENAME_MXE 	= ".MXE.txt";
	private static final String FILENAME_A5SS 	= ".A5SS.txt";
	private static final String FILENAME_A3SS 	= ".A3SS.txt";
	private static final String FILENAME_RI 	= ".RI.txt";
		
	/**
	 * @throws IOException 
	 * 
	 */
	public IndexCreator(String basedir, boolean alsoNovel, String outputFile) throws IOException {
		// try to find files
		File base = new File(basedir);
		if(base.exists()) {
			// open output file
			File file = new File(outputFile);
			FileWriter fw = new FileWriter(file.getAbsoluteFile());
			BufferedWriter bw = new BufferedWriter(fw);
			
			
			ArrayList<MatsOutputParser> parser = new ArrayList<MatsOutputParser>();
			
			// add parsers
			File file_se = new File(base.getAbsolutePath() + File.separator + BASE_NAME + FILENAME_SE);
			File file_mxe = new File(base.getAbsolutePath() + File.separator + BASE_NAME + FILENAME_MXE);
			File file_a5ss = new File(base.getAbsolutePath() + File.separator + BASE_NAME + FILENAME_A5SS);
			File file_a3ss = new File(base.getAbsolutePath() + File.separator + BASE_NAME + FILENAME_A3SS);
			File file_ri = new File(base.getAbsolutePath() + File.separator + BASE_NAME + FILENAME_RI);
			
			if(file_se.exists())
				parser.add(new SEParser(file_se.getAbsolutePath(), bw, null));
			if(file_mxe.exists())
				parser.add(new MXEParser(file_mxe.getAbsolutePath(), bw, null));
			if(file_a5ss.exists())
				parser.add(new A5SSParser(file_a5ss.getAbsolutePath(), bw, null));
			if(file_a3ss.exists())
				parser.add(new A3SSParser(file_a3ss.getAbsolutePath(), bw, null));
			if(file_ri.exists())
				parser.add(new RIParser(file_ri.getAbsolutePath(), bw, null));
			
			// add novel parsers if wished
			if(alsoNovel) {
				File file_se_novel = new File(base.getAbsolutePath() + File.separator + BASE_NAME + NOVE_NAME + FILENAME_SE);
				File file_mxe_novel = new File(base.getAbsolutePath() + File.separator + BASE_NAME + NOVE_NAME + FILENAME_MXE);
				File file_a5ss_novel = new File(base.getAbsolutePath() + File.separator + BASE_NAME + NOVE_NAME + FILENAME_A5SS);
				File file_a3ss_novel = new File(base.getAbsolutePath() + File.separator + BASE_NAME + NOVE_NAME + FILENAME_A3SS);
				File file_ri_novel = new File(base.getAbsolutePath() + File.separator + BASE_NAME + NOVE_NAME + FILENAME_RI);
				
				if(file_se_novel.exists())
					parser.add(new SEParser(file_se_novel.getAbsolutePath(), bw, null));
				if(file_mxe_novel.exists())
					parser.add(new MXEParser(file_mxe_novel.getAbsolutePath(), bw, null));
				if(file_a5ss_novel.exists())
					parser.add(new A5SSParser(file_a5ss_novel.getAbsolutePath(), bw, null));
				if(file_a3ss_novel.exists())
					parser.add(new A3SSParser(file_a3ss_novel.getAbsolutePath(), bw, null));
				if(file_ri_novel.exists())
					parser.add(new RIParser(file_ri_novel.getAbsolutePath(), bw, null));
			}
			// close file handle
			bw.flush();
			bw.close();
			
			// check, if some parsers were found 
			if(parser.size() == 0) {
				file.delete();
			}
		}
	}
	
	/*public static void main(String[] args) throws IOException {
		IndexCreator ic = new IndexCreator("/storageNGS/ngs1/projects/other/Mouse_beckers_huypens/michael/mRNA/MATS/ASEvents", false, "/tmp/fuuu");
	}*/
}
