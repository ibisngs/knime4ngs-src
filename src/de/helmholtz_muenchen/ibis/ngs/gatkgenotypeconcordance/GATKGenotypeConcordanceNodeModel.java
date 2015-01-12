package de.helmholtz_muenchen.ibis.ngs.gatkgenotypeconcordance;

import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.workflow.FlowVariable;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode.GATKNodeModel;

/**
 * This is the model implementation of GATKGenotypeConcordance.
 * 
 *
 * @author Maximilan Hastreiter
 */
public class GATKGenotypeConcordanceNodeModel extends GATKNodeModel {

	
	public static final String CFGKEY_EVAL 		= "EVAL";
	public static final String CFGKEY_COMP 		= "COMP";
	public static final String CFGKEY_UNSAFE	= "UNSAFE";
	
    private final SettingsModelString m_EVAL 	= new SettingsModelString(GATKGenotypeConcordanceNodeModel.CFGKEY_EVAL, "");
    private final SettingsModelString m_COMP 	= new SettingsModelString(GATKGenotypeConcordanceNodeModel.CFGKEY_COMP, "");
	private final SettingsModelBoolean m_UNSAFE = new SettingsModelBoolean(GATKGenotypeConcordanceNodeModel.CFGKEY_UNSAFE,false);
    
	protected GATKGenotypeConcordanceNodeModel(int INPORTS, int OUTPORTS) {
		super(INPORTS, OUTPORTS);
	}

	@Override
	protected String getCommandParameters(final BufferedDataTable[] inData) {
		String command = "--eval "+m_EVAL.getStringValue();
		command 	  += " --comp "+m_COMP.getStringValue();
	
		if(m_UNSAFE.getBooleanValue()){
			command		  += " --unsafe ALL";	
		}
		
		
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
		return m_EVAL.getStringValue()+"_evaluation";
	}

	
	
	@Override
	protected void saveExtraSettingsTo(NodeSettingsWO settings) {
		m_COMP.saveSettingsTo(settings);
		m_EVAL.saveSettingsTo(settings);
		m_UNSAFE.saveSettingsTo(settings);
		
	}

	@Override
	protected void loadExtraValidatedSettingsFrom(NodeSettingsRO settings) throws InvalidSettingsException {
		m_COMP.loadSettingsFrom(settings);
		m_EVAL.loadSettingsFrom(settings);
		m_UNSAFE.loadSettingsFrom(settings);
		
	}

	@Override
	protected void validateExtraSettings(NodeSettingsRO settings) throws InvalidSettingsException {
		m_COMP.validateSettings(settings);
		m_EVAL.validateSettings(settings);
		m_UNSAFE.validateSettings(settings);
		
	}


}

