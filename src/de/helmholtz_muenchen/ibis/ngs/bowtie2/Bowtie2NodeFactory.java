package de.helmholtz_muenchen.ibis.ngs.bowtie2;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTENodeView;

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
        return new HTENodeView<Bowtie2NodeModel>(nodeModel);
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

