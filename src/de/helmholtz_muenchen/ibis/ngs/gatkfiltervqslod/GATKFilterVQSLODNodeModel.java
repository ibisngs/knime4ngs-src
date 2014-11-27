package de.helmholtz_muenchen.ibis.ngs.gatkfiltervqslod;


import de.helmholtz_muenchen.ibis.utils.abstractNodes.SelectVariants.SelectVariantsNodeModel;

/**
 * This is the model implementation of GATKFilterVQSLOD.
 * 
 *
 * @author Maximilian Hastreiter
 */
public class GATKFilterVQSLODNodeModel extends SelectVariantsNodeModel {

	protected GATKFilterVQSLODNodeModel(int INPORTS, int OUTPORTS) {
		super(INPORTS, OUTPORTS);
	
	}

	@Override
	protected String getCommandParameters() {
		
		double Cutoff = Double.parseDouble(getFILTERSTRINGModel().getStringValue());
		
		return "--select_expressions VQSLOD>="+Cutoff;
	}

	@Override
	protected String getOutfileSuffix() {
		return "VQSLOD_Filtered.vcf";
	}
}

