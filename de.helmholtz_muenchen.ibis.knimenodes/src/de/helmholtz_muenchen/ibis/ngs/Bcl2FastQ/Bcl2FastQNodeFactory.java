package de.helmholtz_muenchen.ibis.ngs.Bcl2FastQ;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "Bcl2FastQ" Node.
 * 
 *
 * @author Kaarin Ahomaa
 */
public class Bcl2FastQNodeFactory 
        extends NodeFactory<Bcl2FastQNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public Bcl2FastQNodeModel createNodeModel() {
        return new Bcl2FastQNodeModel();
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
    public NodeView<Bcl2FastQNodeModel> createNodeView(final int viewIndex,
            final Bcl2FastQNodeModel nodeModel) {
        return new Bcl2FastQNodeView(nodeModel);
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
        return new Bcl2FastQNodeDialog();
    }

}

