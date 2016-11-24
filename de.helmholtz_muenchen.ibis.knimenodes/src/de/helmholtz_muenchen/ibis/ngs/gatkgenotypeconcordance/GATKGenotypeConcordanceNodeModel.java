/**
 *  Copyright (C) 2016 the Knime4NGS contributors.
 *  Website: http://ibisngs.github.io/knime4ngs
 *  
 *  This file is part of the KNIME4NGS KNIME extension.
 *  
 *  The KNIME4NGS extension is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


package de.helmholtz_muenchen.ibis.ngs.gatkgenotypeconcordance;


import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataType;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.InvalidSettingsException;

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

