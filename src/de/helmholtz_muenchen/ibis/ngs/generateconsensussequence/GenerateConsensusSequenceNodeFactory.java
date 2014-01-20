package de.helmholtz_muenchen.ibis.ngs.generateconsensussequence;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "GenerateConsensusSequence" Node.
 * 
 *
 * @author 
 */
public class GenerateConsensusSequenceNodeFactory 
        extends NodeFactory<GenerateConsensusSequenceNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public GenerateConsensusSequenceNodeModel createNodeModel() {
        return new GenerateConsensusSequenceNodeModel();
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
    public NodeView<GenerateConsensusSequenceNodeModel> createNodeView(final int viewIndex,
            final GenerateConsensusSequenceNodeModel nodeModel) {
        return new GenerateConsensusSequenceNodeView(nodeModel);
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
        return new GenerateConsensusSequenceNodeDialog();
    }

}

