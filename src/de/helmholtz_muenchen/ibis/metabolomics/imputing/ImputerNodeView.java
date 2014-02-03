package de.helmholtz_muenchen.ibis.metabolomics.imputing;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RNodeView;

/**
 * <code>NodeView</code> for the "Imputer" Node.
 * 
 *
 * @author Jonas Zierer
 */
public class ImputerNodeView extends RNodeView<ImputerNodeModel> {

    /**
     * Creates a new view.
     * 
     * @param nodeModel The model (class: {@link ImputerNodeModel})
     */
    protected ImputerNodeView(final ImputerNodeModel nodeModel) {
        super(nodeModel);
    }

}

