package de.helmholtz_muenchen.ibis.ngs.gatkexcludevariants;

import org.knime.core.node.InvalidSettingsException;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.SelectVariants.SelectVariantsNodeModel;

/**
 * This is the model implementation of GATKExcludeVariants.
 * 
 *
 * @author Tim Jeske
 */
public class GATKExcludeVariantsNodeModel extends SelectVariantsNodeModel {
    
    /**
     * Constructor for the node model.
     */
    protected GATKExcludeVariantsNodeModel() {
        super(1, 1);
        getFILTERSTRINGModel().setStringValue("--excludeFiltered --excludeNonVariants");
    }

	@Override
	protected String getCommandParameters() {
		return getFILTERSTRINGModel().getStringValue();
	}

	@Override
	protected String getOutfileSuffix() {
		return ".PASS_variants.vcf";
	}

	@Override
	protected void extraConfig() throws InvalidSettingsException {
		// TODO Auto-generated method stub
		
	}
}

