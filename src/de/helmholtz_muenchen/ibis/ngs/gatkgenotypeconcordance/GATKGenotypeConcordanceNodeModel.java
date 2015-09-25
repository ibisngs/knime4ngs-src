package de.helmholtz_muenchen.ibis.ngs.gatkgenotypeconcordance;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;

import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode.GATKNodeModel;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;

/**
 * This is the model implementation of GATKGenotypeConcordance.
 * 
 *
 * @author Maximilan Hastreiter
 */
public class GATKGenotypeConcordanceNodeModel extends GATKNodeModel {

	
	public static final String CFGKEY_EVAL 		= "EVAL";
	public static final String CFGKEY_COMP 		= "COMP";
//	public static final String CFGKEY_UNSAFE	= "UNSAFE";
	
    private final SettingsModelString m_EVAL 	= new SettingsModelString(GATKGenotypeConcordanceNodeModel.CFGKEY_EVAL, "");
    private final SettingsModelString m_COMP 	= new SettingsModelString(GATKGenotypeConcordanceNodeModel.CFGKEY_COMP, "");
//	private final SettingsModelBoolean m_UNSAFE = new SettingsModelBoolean(GATKGenotypeConcordanceNodeModel.CFGKEY_UNSAFE,false);
    
    private String OUTFILE, LOCKFILE;
    
	protected GATKGenotypeConcordanceNodeModel() {
		super(OptionalPorts.createOPOs(2,1,2), OptionalPorts.createOPOs(1));
	}

	@Override
	protected String getCommandParameters(final BufferedDataTable[] inData) throws InvalidSettingsException {
		
		String eval, comp;
		
		if (inData[0]!=null) {
			eval = inData[0].iterator().next().getCell(0).toString();
		} else {
			eval = m_EVAL.getStringValue();
		}
		
		if (inData[1]!=null) {
			comp = inData[1].iterator().next().getCell(0).toString();
		} else {
			comp = m_COMP.getStringValue();
		}
		
		if((eval.equals("") || Files.notExists(Paths.get(eval))) || (comp.equals("") || Files.notExists(Paths.get(comp)))) {
			throw new InvalidSettingsException("Specify two VCF files as a test and a truth set!");
		}
		
		OUTFILE = eval + "_evaluation";
		LOCKFILE = IO.replaceFileExtension(eval, SuccessfulRunChecker.LOCK_ENDING);
		
		String command = "--eval "+eval;
		command 	  += " --comp "+comp;
	
//		if(m_UNSAFE.getBooleanValue()){
//			command		  += " --unsafe ALL";	
//		}
		
		
		//Push FlowVars in order to provide Infile Names for plotgenotypeconcordance Node
		pushFlowVariableString("EVAL", m_EVAL.getStringValue());
		pushFlowVariableString("COMP", m_COMP.getStringValue());
		
		return command;
	}

	@Override
	protected String getCommandWalker() {
		return "GenotypeConcordance";
	}

	@Override
	protected String getOutfile() {
		return OUTFILE;
	}

	
	
	@Override
	protected void saveExtraSettingsTo(NodeSettingsWO settings) {
		m_COMP.saveSettingsTo(settings);
		m_EVAL.saveSettingsTo(settings);
//		m_UNSAFE.saveSettingsTo(settings);
		
	}

	@Override
	protected void loadExtraValidatedSettingsFrom(NodeSettingsRO settings) throws InvalidSettingsException {
		m_COMP.loadSettingsFrom(settings);
		m_EVAL.loadSettingsFrom(settings);
//		m_UNSAFE.loadSettingsFrom(settings);
		
	}

	@Override
	protected void validateExtraSettings(NodeSettingsRO settings) throws InvalidSettingsException {
		m_COMP.validateSettings(settings);
		m_EVAL.validateSettings(settings);
//		m_UNSAFE.validateSettings(settings);
		
	}

	@Override
	protected File getLockFile() {
		return new File(LOCKFILE);
	}
}

