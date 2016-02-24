package de.helmholtz_muenchen.ibis.ngs.gatkgenotypeconcordance;


import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataType;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode.GATKNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;

/**
 * This is the model implementation of GATKGenotypeConcordance.
 * 
 *
 * @author Maximilan Hastreiter
 */
public class GATKGenotypeConcordanceNodeModel extends GATKNodeModel {

    private String OUTFILE;
    private int eval_index, comp_index;
    
	protected GATKGenotypeConcordanceNodeModel() {
		super(OptionalPorts.createOPOs(2), OptionalPorts.createOPOs(1));
	}

	@Override
	protected String getCommandParameters(final BufferedDataTable[] inData) throws InvalidSettingsException {
		
		String eval, comp;
		
    	try{
    		eval = inData[0].iterator().next().getCell(eval_index).toString();
    		if(!eval.endsWith(".vcf")){
    			throw new InvalidSettingsException("A cell of the first input table has to be the path to the evaluation VCF infile but it is "+eval);
    		}
    	}catch(IndexOutOfBoundsException e){
    			throw new InvalidSettingsException("A cell of the first input table has to be the path to the evaluation VCF infile but it is empty.");
    	}
		
    	try{
    		comp = inData[1].iterator().next().getCell(comp_index).toString();
    		if(!comp.endsWith(".vcf")){
    			throw new InvalidSettingsException("A cell of the second input table has to be the path to the comparison VCF infile but it is "+comp);
    		}
    	}catch(IndexOutOfBoundsException e){
    			throw new InvalidSettingsException("A cell of the second input table has to be the path to the comparison VCF infile but it is empty.");
    	}
		
		
		OUTFILE = eval + "_evaluation";
		
		
		String command = "--eval "+eval;
		command 	  += " --comp "+comp;
		
		//Push FlowVars in order to provide Infile Names for plotgenotypeconcordance Node
		pushFlowVariableString("EVAL", eval);
		pushFlowVariableString("COMP", comp);
		
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
	
	}

	@Override
	protected void loadExtraValidatedSettingsFrom(NodeSettingsRO settings) throws InvalidSettingsException {
	
	}

	@Override
	protected void validateExtraSettings(NodeSettingsRO settings) throws InvalidSettingsException {
	
	}


	@Override
	protected boolean checkInputCellType(DataTableSpec[] inSpecs) {
		comp_index = -1;
		eval_index = -1;
		
		for(int i = 0; i < inSpecs[0].getNumColumns(); i++) {
    		if(inSpecs[0].getColumnSpec(i).getType().toString().equals("VCFCell")) {
    			eval_index = i;
    		}
    	}
		
		for(int i = 0; i < inSpecs[1].getNumColumns(); i++) {
    		if(inSpecs[1].getColumnSpec(i).getType().toString().equals("VCFCell")) {
    			comp_index = i;
    		}
    	}
		return (comp_index>-1 && eval_index>-1);
	}

	@Override
	protected DataType getOutColType() {
		return FileCell.TYPE;
	}

	@Override
	protected void extraConfig() throws InvalidSettingsException {
		// TODO Auto-generated method stub
		
	}
}

