package de.helmholtz_muenchen.ibis.ngs.snpeff;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "SnpEff" Node.
 * 
 *
 * @author 
 */
public class SnpEffNodeFactory 
        extends NodeFactory<SnpEffNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public SnpEffNodeModel createNodeModel() {
        return new SnpEffNodeModel();
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
    public NodeView<SnpEffNodeModel> createNodeView(final int viewIndex,
            final SnpEffNodeModel nodeModel) {
        return new SnpEffNodeView(nodeModel);
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
        return new SnpEffNodeDialog();
    }

}

