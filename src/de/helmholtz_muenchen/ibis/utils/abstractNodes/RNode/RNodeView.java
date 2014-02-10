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
package de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode;

import java.awt.Dimension;

import org.knime.core.node.NodeView;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.ScriptNode.LogPanel;

/**
 * <code>NodeView</code> for the "RNodeModel" Node.
 * 
 * @author Jonas Zierer
 */
public abstract class RNodeView<T extends RNodeModel> extends NodeView<T> {
    

    // panel which actually paints the bins
	private final LogPanel m_panel_log;
    
    /**
     * Creates a new view.
     * 
     * @param nodeModel the model class: {@link RNodeModel}
     */
	protected RNodeView(final T nodeModel) {
        super(nodeModel);
        m_panel_log = new LogPanel(nodeModel.getSTDOUT(), nodeModel.getSTDERR());
        m_panel_log.setPreferredSize(new Dimension(800, 600));
        setComponent(m_panel_log);    
        this.setViewTitleSuffix("log");
        
    }

    /** {@inheritDoc} */
    @Override
    protected void modelChanged() {
    	m_panel_log.updateView();
    }

    /** {@inheritDoc} */
    @Override
    protected void onClose() {

    }

    /** {@inheritDoc} */
    @Override
    protected void onOpen() {
 
    }



}
