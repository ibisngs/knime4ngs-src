package de.helmholtz_muenchen.ibis.ngs.impute2;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "Impute2" Node.
 * 
 *
 * @author Tanzeem Haque
 */
public class Impute2NodeFactory 
        extends NodeFactory<Impute2NodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public Impute2NodeModel createNodeModel() {
        return new Impute2NodeModel();
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
    public NodeView<Impute2NodeModel> createNodeView(final int viewIndex,
            final Impute2NodeModel nodeModel) {
        return new Impute2NodeView(nodeModel);
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
        return new Impute2NodeDialog();
    }

}

