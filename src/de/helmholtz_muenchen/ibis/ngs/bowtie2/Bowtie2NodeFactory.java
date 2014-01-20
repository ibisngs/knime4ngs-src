package de.helmholtz_muenchen.ibis.ngs.bowtie2;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "Bowtie2" Node.
 * 
 *
 * @author 
 */
public class Bowtie2NodeFactory 
        extends NodeFactory<Bowtie2NodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public Bowtie2NodeModel createNodeModel() {
        return new Bowtie2NodeModel();
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
    public NodeView<Bowtie2NodeModel> createNodeView(final int viewIndex,
            final Bowtie2NodeModel nodeModel) {
        return new Bowtie2NodeView(nodeModel);
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
        return new Bowtie2NodeDialog();
    }

}

