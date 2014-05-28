package de.helmholtz_muenchen.ibis.misc.regression;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RNodeView;

/**
 * <code>NodeFactory</code> for the "OutlierRemoval" Node.
 * Remove Outliers from Data Columns
 *
 * @author Jonas Zierer
 */
public class RegressionNodeFactory  extends NodeFactory<RegressionNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public RegressionNodeModel createNodeModel() {
        return new RegressionNodeModel();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getNrNodeViews() {
        return 1;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public NodeView<RegressionNodeModel> createNodeView(final int viewIndex, final RegressionNodeModel nodeModel) {
        return new RNodeView<RegressionNodeModel>(nodeModel);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean hasDialog() {
        return true;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public NodeDialogPane createNodeDialogPane() {
        return new RegressionNodeDialog();
    }

}

