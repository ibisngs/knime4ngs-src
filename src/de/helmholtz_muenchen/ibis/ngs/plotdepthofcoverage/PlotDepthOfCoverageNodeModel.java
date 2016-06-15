package de.helmholtz_muenchen.ibis.ngs.plotdepthofcoverage;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;

import org.knime.core.data.DataRow;
import org.knime.core.data.DataTableSpec;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RNodeModel;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;

/**
 * This is the model implementation of PlotDepthOfCoverage.
 * 
 *
 * @author Tim Jeske
 */
public class PlotDepthOfCoverageNodeModel extends RNodeModel {
    
	protected static final NodeLogger LOGGER = NodeLogger.getLogger(PlotDepthOfCoverageNodeModel.class);
	
	/**
	 * Config Keys
	 */
	
	public static final String CFGKEY_INFOLDER = "infolder";
	public static final String CFGKEY_FILESUFFIX = "suffix";
	
	/**
	 * Node Models
	 */
//	private final SettingsModelString m_infolder = new SettingsModelString(PlotDepthOfCoverageNodeModel.CFGKEY_INFOLDER,"");
//	private final SettingsModelString m_filesuffix = new SettingsModelString(PlotDepthOfCoverageNodeModel.CFGKEY_FILESUFFIX,"");
	
	private static final String SCRIPT_PATH = "ngs" + File.separatorChar + "de" + File.separatorChar + "plotExomeCapture.R";
	
//	public boolean optionalPort=false;
	
    /**
     * Constructor for the node model.
     */
    protected PlotDepthOfCoverageNodeModel() {
    	super(OptionalPorts.createOPOs(1), OptionalPorts.createOPOs(0), SCRIPT_PATH, new String[]{"--infile"}, new String[]{});
//    	addSetting(m_filesuffix);
//    	addSetting(m_infolder);
    }

    
    /**
     * {@inheritDoc}
     * @throws Exception 
     */
    @Override
	protected BufferedDataTable[] execute(final BufferedDataTable[] inData, final ExecutionContext exec) throws Exception{    	
    	
    	ArrayList<String> files = new ArrayList<>();
    	String outfolder;
    	
//    	if(optionalPort) {
    		String firstIn = inData[0].iterator().next().getCell(0).toString();
    		outfolder = firstIn.substring(0,firstIn.lastIndexOf(File.separatorChar));
    		Iterator<DataRow> iter = inData[0].iterator();
    		while(iter.hasNext()) {
    			files.add(iter.next().getCell(0).toString());
    		}
//    	} else {
//    		outfolder = m_infolder.getStringValue();
//    		FileSearch fileSearch = new FileSearch(); 
//
//        	//try different directory and filename :)
//        	fileSearch.searchDirectory(new File(m_infolder.getStringValue()), m_filesuffix.getStringValue());
//    	    for (String matched : fileSearch.getResult()){
//    	    	LOGGER.debug("Found : " + matched);
//    	    	files.add(matched);
//    	    }
//    	}

    	this.addArgument("--in"   , writeDocSampleSummary(outfolder,files));
    	super.execute(inData, exec);
    	
		return null;
	}
    
    public String writeDocSampleSummary (String outfolder, ArrayList<String> infiles) throws IOException {
    	String outfile;
    	if(outfolder.endsWith(File.separatorChar+"")) {
    		outfile = outfolder + "DoCSampleSummary.csv";
    	} else {
    		outfile = outfolder + File.separatorChar + "DoCSampleSummary.csv";
    	}
    	
    	PrintWriter writer = new PrintWriter(outfile, "UTF-8");
    	for(String file: infiles) {
	    		String everything ="";
	    		try {
	    			everything = processCoverageFile(file);
	    		} catch (IOException e) {
	    			setWarningMessage("Could not process file: "+file);
	    		}
	    		writer.println(file+"\t"+everything);
    	}
    	writer.close();
    	return outfile;
    }
    
	private String processCoverageFile(String DoCOutfile) throws IOException{

		String result = "";

		int[] Coverage = new int[102];

		BufferedReader br = new BufferedReader(new FileReader(DoCOutfile));

		String line;
		double totalBasePairs = 0;

		br.readLine(); // Skip Header

		while ((line = br.readLine()) != null) {
			totalBasePairs++;
			int curr_coverage = Integer.parseInt(line.split("\t")[3]);

			if (curr_coverage > 101) {
				curr_coverage = 101;
			}

			Coverage[curr_coverage]++;
		}

		double remaining = totalBasePairs;
		result = "100.00";
		for (int i = 1; i < 101; i++) {
			remaining = remaining - Coverage[i - 1]; // Anzahl von BasePairs mit
														// Cov > Coverage[i]
			double value = remaining / totalBasePairs * 100;
			value = (double) Math.round(value * 100) / 100;
			result += "\t" + value;
		}
		br.close();
		return result;
	}

    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {

    	return null;
    }
}

