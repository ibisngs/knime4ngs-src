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

package de.helmholtz_muenchen.ibis.ngs.vat;


import org.knime.core.node.ExecutionContext;

import de.helmholtz_muenchen.ibis.utils.threads.Executor;

/**
 * 
 * @author Marie-Sophie Friedl
 */

public class RunVAT {
	
	public static void IndelMapper(ExecutionContext exec, String indelmapper, String intervals, String transcripts, String vcfin, String vcfout) throws Exception{
		
		String cmd=indelmapper;
		cmd+=" "+intervals;
		cmd+=" "+transcripts;
		
		VATNodeModel.logger.info("Running VAT indelMapper...");
		VATNodeModel.logger.info("Log file can be found in "+vcfout+".err.log");
		
		Executor.executeCommand(new String []{cmd}, exec, new String[]{"LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/share/tmp/knime/libs/libbios-1.0.0/lib:/home/share/tmp/knime/libs/gsl-1.16/lib:/home/share/tmp/knime/libs/gd/2.0.35/lib"}, VATNodeModel.logger, vcfout, vcfout+".err.log", vcfin);
		
	}
	
	public static void SNPMapper(ExecutionContext exec, String snpmapper, String intervals, String transcripts, String vcfin, String vcfout)throws Exception{
		
		String cmd=snpmapper;
		cmd+=" "+intervals;
		cmd+=" "+transcripts;
		
		VATNodeModel.logger.info("Running VAT snpMapper...");
		VATNodeModel.logger.info("Log file can be found in "+vcfout+".err.log");
		
		Executor.executeCommand(new String []{cmd}, exec, new String[]{"LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/share/tmp/knime/libs/libbios-1.0.0/lib:/home/share/tmp/knime/libs/gsl-1.16/lib:/home/share/tmp/knime/libs/gd/2.0.35/lib"}, VATNodeModel.logger, vcfout, vcfout+".err.log", vcfin);

	}
	
}