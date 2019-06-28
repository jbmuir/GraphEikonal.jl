module FEFMM_Tests
    using Test
    using GraphEikonal
    using LightGraphs
    using SimpleWeightedGraphs

    a23 = (0.8/0.9)^2
    a45 = ((1.5-0.9*sqrt(0.01)-1.8*sqrt(0.02)-3*sqrt(0.03)-2.2*sqrt(0.04))/1.3)^2
    a56 = (0.7/0.9)^2
    A = [1.0 0.16 0.0 0.03 0.0 6.25;
         0.16 1.0 a23 0.02 0.0 1.1025;
         0.0 a23 1.0 0.01 0.0 0.0;
         0.03 0.02 0.01 1.0 a45 0.04;
         0.0 0.0 0.0 a45 1.0 a56;
         6.25 1.1025 0.0 0.04 a56 1.0]

    b23 = (0.8/0.9)
    b26 = sqrt(0.9^2-0.4^2*1.2^2)/0.4
    b45 = sqrt(1.5^2-0.02^2*3^2-0.01^2*1.8^2-0.03^2*0.9^2-0.04^2*2.2^2)/1.3
    b56 = (0.7/0.9)
    B = [1.0 0.4 0.0 0.02 0.0 2.5;
          0.4 1.0 b23 0.01 0.0 b26;
          0.0 b23 1.0 0.03 0.0 0.0;
          0.02 0.01 0.03 1.0 b45 0.04;
          0.0 0.0 0.0 b45 1.0 b56;
          2.5 b26 0.0 0.04 b56 1.0]

    c23 = (0.8/0.9)^2
    c34 = (1.5/0.9)^2
    c45 = (0.7/0.9)^2
    C = [1.0 0.5625 0.0 0.03 0.0 6.25;
         0.5625 1.0 c23 0.02 0.0 0.01;
         0.0 c23 1.0 c34 0.0 0.0;
         0.03 0.02 c34 1.0 c45 0.04;
         0.0 0.0 0.0 c45 1.0 a56;
         6.25 0.01 0.0 0.04 a56 1.0]


    p = [1.1, 0.9, 0.8, 1.5, 0.7, 2.0]
    sol = [0.0, 1.2, 2.1, 3.0, 1.7, 0.8]
    g1 = SimpleWeightedGraph(A)
    g2 = SimpleWeightedGraph(B.^2)
    g3 = SimpleWeightedGraph(C)

    @testset "Graph Tests" begin
    @test all(isapprox.(graph_eikonal(g1,[1]; slowv=p, norm=L1Norm()), sol))
    @test all(isapprox.(graph_eikonal(g2,[1]; slowv=p, norm=L2Norm()), sol))
    @test all(isapprox.(graph_eikonal(g3,[1]; slowv=p, norm=LInfNorm()), sol))

    end
    
end