package de.helmholtz_muenchen.ibis.ngs.gatkselectvariants;


import org.knime.core.node.InvalidSettingsException;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.SelectVariants.SelectVariantsNodeModel;



/**
 * This is the model implementation of GATKSelectVariants.
 * 
 *
 * @author Maximilian Hastreiter
 */
public class GATKSelectVariantsNodeModel extends SelectVariantsNodeModel {

	
    protected GATKSelectVariantsNodeModel(int INPORTS, int OUTPORTS) {
		super(INPORTS, OUTPORTS);
	}

	@Override
	protected String getCommandParameters() {
		return "-selectType "+getFILTERSTRINGModel().getStringValue();
	}

	@Override
	protected String getOutfileSuffix() {
		return "."+getFILTERSTRINGModel().getStringValue()+".vcf";
	}

	@Override
	protected void extraConfig() throws InvalidSettingsException {
		// TODO Auto-generated method stub
		
	}    
}


