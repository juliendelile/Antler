import sledge

sledge.generate_dataset(
    n_cells                  = 2000,
    n_tree_states            = [19, 4],
    p_branching              = [0.4, 0],
    n_genes_per_tree_state   = [80, 80],
    n_cc_states              = 4,
    n_genes_per_cc_state     = 100,
    n_unexpressed_genes      = 1000,
    noise_intensity          = 0.3,
    seed                     = 13579,
    plot_size_factor         = 1)


