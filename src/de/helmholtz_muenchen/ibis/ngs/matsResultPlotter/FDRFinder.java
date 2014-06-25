package de.helmholtz_muenchen.ibis.ngs.matsResultPlotter;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import de.helmholtz_muenchen.ibis.ngs.matsResultIndexer.A3SSParser;
import de.helmholtz_muenchen.ibis.ngs.matsResultIndexer.A5SSParser;
import de.helmholtz_muenchen.ibis.ngs.matsResultIndexer.MXEParser;
import de.helmholtz_muenchen.ibis.ngs.matsResultIndexer.MatsOutputParser;
import de.helmholtz_muenchen.ibis.ngs.matsResultIndexer.RIParser;
import de.helmholtz_muenchen.ibis.ngs.matsResultIndexer.SEParser;

/**
 * Finds 
 * @author Michael Kluge
 *
 */
public class FDRFinder {
	
	private static final String JUNCTION_NAME	= "JunctionCountOnly.txt";
	private static final String READS_NAME 		= "ReadsOnTargetAndJunctionCounts.txt";
	
	private static final String FILENAME_SE 	= "SE.MATS.";
	private static final String FILENAME_MXE 	= "MXE.MATS.";
	private static final String FILENAME_A5SS 	= "A5SS.MATS.";
	private static final String FILENAME_A3SS 	= "A3SS.MATS.";
	private static final String FILENAME_RI 	= "RI.MATS.";
	
	private final HashMap<String, ArrayList<String>> EVENTS = new HashMap<String, ArrayList<String>>(); 
		
	/**
	 * @throws IOException 
	 * 
	 */
	public FDRFinder(String basedir, boolean alsoReads, double FDR) throws IOException {
		// try to find files
		File base = new File(basedir);
		if(base.exists()) {
			ArrayList<MatsOutputParser> parser = new ArrayList<MatsOutputParser>();
			
			// add only junction files
			if(!alsoReads) {
				// add parsers
				File file_se = new File(base.getAbsolutePath() + File.separator + FILENAME_SE + JUNCTION_NAME);
				File file_mxe = new File(base.getAbsolutePath() + File.separator + FILENAME_MXE + JUNCTION_NAME);
				File file_a5ss = new File(base.getAbsolutePath() + File.separator + FILENAME_A5SS + JUNCTION_NAME);
				File file_a3ss = new File(base.getAbsolutePath() + File.separator + FILENAME_A3SS + JUNCTION_NAME);
				File file_ri = new File(base.getAbsolutePath() + File.separator + FILENAME_RI + JUNCTION_NAME);
				
				if(file_se.exists())
					parser.add(new SEParser(file_se.getAbsolutePath(), null, FDR));
				if(file_mxe.exists())
					parser.add(new MXEParser(file_mxe.getAbsolutePath(), null, FDR));
				if(file_a5ss.exists())
					parser.add(new A5SSParser(file_a5ss.getAbsolutePath(), null, FDR));
				if(file_a3ss.exists())
					parser.add(new A3SSParser(file_a3ss.getAbsolutePath(), null, FDR));
				if(file_ri.exists())
					parser.add(new RIParser(file_ri.getAbsolutePath(), null, FDR));
			}
			// add read parsers
			else {
				File file_se = new File(base.getAbsolutePath() + File.separator + FILENAME_SE + READS_NAME);
				File file_mxe = new File(base.getAbsolutePath() + File.separator + FILENAME_MXE + READS_NAME);
				File file_a5ss = new File(base.getAbsolutePath() + File.separator + FILENAME_A5SS + READS_NAME);
				File file_a3ss = new File(base.getAbsolutePath() + File.separator + FILENAME_A3SS + READS_NAME);
				File file_ri = new File(base.getAbsolutePath() + File.separator + FILENAME_RI + READS_NAME);
				
				if(file_se.exists())
					parser.add(new SEParser(file_se.getAbsolutePath(), null, FDR));
				if(file_mxe.exists())
					parser.add(new MXEParser(file_mxe.getAbsolutePath(), null, FDR));
				if(file_a5ss.exists())
					parser.add(new A5SSParser(file_a5ss.getAbsolutePath(), null, FDR));
				if(file_a3ss.exists())
					parser.add(new A3SSParser(file_a3ss.getAbsolutePath(), null, FDR));
				if(file_ri.exists())
					parser.add(new RIParser(file_ri.getAbsolutePath(), null, FDR));
			}
			
			// get names
			for(Iterator<MatsOutputParser> it = parser.iterator(); it.hasNext(); ) {
				MatsOutputParser p = it.next();
				this.EVENTS.put(p.getType(), p.getSigEvents());
			}
		}
	}
	
	/**
	 * Returns the significant events for each type
	 * @return
	 */
	public HashMap<String, ArrayList<String>> getSigEvents() {
		return this.EVENTS;
	}
}
