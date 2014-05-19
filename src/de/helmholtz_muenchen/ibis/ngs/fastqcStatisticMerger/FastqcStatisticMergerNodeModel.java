package de.helmholtz_muenchen.ibis.ngs.fastqcStatisticMerger;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.StatisticMerger.StatisticMergerNodeModel;

/**
 * This is the model implementation of FastqcStatisticMergerNodeModel.
 * This node can be used to merge the statistics of the FASTQC Node
 *
 * @author Michael Kluge
 */
public class FastqcStatisticMergerNodeModel extends StatisticMergerNodeModel {
    
	public static final String DATA_FILE = "fastqc_data.txt";
	public static final String MODULE_PREFIX = ">>";
	public static final String MODULE_END = "END_MODULE";
	public static final String BASIC_FILENAME = "Filename";
	public static final String BASIC_MODULE = "Basic Statistics";
	public static final String TAB = "\t";
	public static final String HEADER_PREFIX = "#";
	public static final String HEADER_NAME = "Filename";
	public static final String HEADER_STATUS = "#Module" + TAB + "Status";
	public static final String FILE_ENDING = ".txt";
	protected static final String KEY = "MODULE_KEY_";
	protected static final String MERGER_NAME = "FastQC";
	
	// Available module names
	private static final String MODULE_BASIC = "Basic Statistics";
	private static final String MODULE_BSQ = "Per base sequence quality";
	private static final String MODULE_SQS = "Per sequence quality scores";
	private static final String MODULE_BSC = "Per base sequence content";
	private static final String MODULE_BGCC = "Per base GC content";
	private static final String MODULE_SGCC = "Per sequence GC content";
	private static final String MODULE_BNC = "Per base N content";
	private static final String MODULE_SLD = "Sequence Length Distribution";
	private static final String MODULE_SDL = "Sequence Duplication Levels";
	private static final String MODULE_OS = "Overrepresented sequences";
	private static final String MODULE_FOC = "Filter-Options Collector";
	private static final String MODULE_STATUS = "Merge Status of all Modules";
	protected static final ArrayList<String> MODULE_NAMES = new ArrayList<String>();
	
	 /**
     * add the used settings
     */
    static {
        // add modules
        MODULE_NAMES.add(MODULE_BASIC);
        MODULE_NAMES.add(MODULE_BSQ);
        MODULE_NAMES.add(MODULE_SQS);
        MODULE_NAMES.add(MODULE_BSC);
        MODULE_NAMES.add(MODULE_BGCC);
        MODULE_NAMES.add(MODULE_SGCC);
        MODULE_NAMES.add(MODULE_BNC);
        MODULE_NAMES.add(MODULE_SLD);
        MODULE_NAMES.add(MODULE_SDL);
        MODULE_NAMES.add(MODULE_OS);
        MODULE_NAMES.add(MODULE_FOC);
        MODULE_NAMES.add(MODULE_STATUS); 
        
        addModuleNames(MERGER_NAME, MODULE_NAMES);
    }
	
	@Override
	protected boolean extractTable(String moduleName, String file, BufferedWriter outfile, boolean writeHeader) {
		boolean isStatusMode = MODULE_STATUS.equals(moduleName);
		File f = new File(file);
		// test if correct file is given
		if(!f.isFile() || !f.canRead() || !f.getName().equals(DATA_FILE))
			return false;
		
		// open file
		try {
			String line;
			BufferedReader r = new BufferedReader(new FileReader(f));

			String name = null;
			boolean hasFoundBasicModule = false;
			boolean hasFoundCorrectModule = false;

			// read lines
			while((line = r.readLine()) != null) {
				// find start of basic module
				if(!hasFoundBasicModule && !hasFoundCorrectModule) {
					if(line.startsWith(MODULE_PREFIX + BASIC_MODULE)) 
						hasFoundBasicModule = true;
				}
				// get filename of fasta file which was analysed with FASTQ
				else if(hasFoundBasicModule && !hasFoundCorrectModule) {
					// find name of file
					if(name == null && line.startsWith(BASIC_FILENAME)) {
						name = line.replace(BASIC_FILENAME, "").replace(TAB, "");
						// go back to start of file (if basic module should be merged)
						r = new BufferedReader(new FileReader(f));
						// check, if status should be merged
						if(isStatusMode) 
							hasFoundCorrectModule = true; // all modules are correct ;)
					}
					else if(line.startsWith(MODULE_PREFIX + moduleName)) 
						hasFoundCorrectModule = true;
				}
				// is in correct module --> write lines
				else if(hasFoundCorrectModule) {
					if(!isStatusMode && line.equals(MODULE_PREFIX + MODULE_END)) 
						break;
					else {
						// check, if starts with char for header line
						if(line.startsWith(HEADER_PREFIX)) {
							if(writeHeader) {
								if(!isStatusMode)
									outfile.write(line + TAB + HEADER_NAME);
								else
									outfile.write(HEADER_STATUS + TAB + HEADER_NAME);
								// do not write another header in case of status merger
								outfile.newLine();
								writeHeader = false;
							}
						}
						else if(!isStatusMode) {
							outfile.write(line + TAB + name);
							outfile.newLine();
						}
						// statusMode
						else if(line.startsWith(MODULE_PREFIX) && !line.equals(MODULE_PREFIX + MODULE_END)) {
							String currentModuleName = line.replaceFirst(MODULE_PREFIX, "").split(TAB)[0];
							outfile.write(currentModuleName);
							outfile.write(TAB);
							outfile.write(line.replaceFirst(MODULE_PREFIX + currentModuleName + TAB, "")); 
							outfile.write(TAB);
							outfile.write(name);
							outfile.newLine();
						}
					}
				}
			}
			outfile.flush();
			// close the file
			r.close();
			return true;
		}
		catch(Exception e) {
			e.printStackTrace();
		}
		return false;
	}

	@Override
	public String getNameOfStatisticFile() {
		return DATA_FILE;
	}

	@Override
	public String getMergerName() {
		return MERGER_NAME;
	}

	@Override
	public void finalize(BufferedWriter outfile) throws IOException {}
}

