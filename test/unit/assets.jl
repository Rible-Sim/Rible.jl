using Test
import Rible as RB

@testset "Asset Paths" begin
    # Tier 1 asset (local)
    @test isfile(RB.assetpath("stars.jpg"))
    @test isfile(RB.assetpath("assembly1.STL"))

    # Tier 2 asset (lazy download via Artifacts)
    # This will trigger the artifact download on first run
    @test isfile(RB.assetpath("LJ.STL"))
    @test isfile(RB.assetpath("Toupise.STL"))
end
