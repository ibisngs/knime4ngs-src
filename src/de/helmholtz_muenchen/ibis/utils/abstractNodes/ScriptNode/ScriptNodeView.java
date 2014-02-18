/*
 * ------------------------------------------------------------------------
  * This source code, its documentation and all appendant files
 * are protected by copyright law. All rights reserved.
 *
 * Copyright, 2008 - 2012
 * KNIME.com, Zurich, Switzerland
 *
 * You may not modify, publish, transmit, transfer or sell, reproduce,
 * create derivative works from, distribute, perform, display, or in
 * any way exploit any of the content, in whole or in part, except as
 * otherwise expressly permitted in writing by the copyright owner or
 * as specified in the license file distributed with this product.
 *
 * If you have any questions please contact the copyright holder:
 * website: www.knime.com
 * email: contact@knime.com
 * ---------------------------------------------------------------------
 */
package de.helmholtz_muenchen.ibis.utils.abstractNodes.ScriptNode;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.ExecutorNode.ExecutorNodeView;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RNodeModel;


/**
 * <code>NodeView</code> for the "RNodeModel" Node.
 * 
 * @author Jonas Zierer
 */
public abstract class ScriptNodeView<T extends ScriptNodeModel> extends ExecutorNodeView<T> {
    
	/**
     * Creates a new view.
     * 
     * @param nodeModel the model class: {@link RNodeModel}
     */
	protected ScriptNodeView(final T nodeModel) {
        super(nodeModel);
    }
}
