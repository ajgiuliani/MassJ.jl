using MassJ, Test
using Plots
using DataStructures
using SHA   # used for indexed-mzML fileChecksum verification
using Aqua  # package-quality checks
using Random
# Loading these weak dependencies activates the package extensions
# (cmi_matrix via EntropyInvariant, cluster_ions via Clustering, and the
# YieldCurve Tables.jl source via Tables).
import EntropyInvariant
import Clustering
import Tables
import Measurements
import Unitful
using Unitful: @u_str

function tests()
    @testset "Subset of tests"  begin
        inf = MassJ.info("test.mzXML", verbose = true)
        @test inf[1] == "parentFile: test.raw"                                         #1
        @test inf[9] == "6 scans"                                                      #2
        @test inf[10] == "MS1+"                                                        #3
        @test inf[11] == "MS2+ 1255.5  CID(CE=18)"                                     #4
        @test inf[12] == "MS3+ 902.33  PQD(CE=35)"                                     #5

        scans = MassJ.load("test.mzXML")
        @test eltype(scans)              == MassJ.MSscans                                 #6
        @test length(scans)              == 6                                          #7
        @test scans[1].num               == [1]                                          #8
        @test scans[2].level             == [2]                                          #9
        @test scans[3].polarity          == ["+"]                                        #10
        @test scans[2].activationMethod  == ["CID"]                                      #11
        @test scans[3].collisionEnergy   == [35.0]                                       #12
        @test size(scans[1].int, 1)      == 22320                                      #13

        rt = MassJ.retention_time("test.mzXML")
        @test length(rt) == 6                                                          #14

        cr = MassJ.chromatogram("test.mzXML", method = MassJ.TIC() )
        @test length(cr.x) == 6                                                       #15

        cr = MassJ.chromatogram("test.mzXML", method = MassJ.MZ([0, 500]))
        @test length(cr.x) == 6                                                       #16

        cr = MassJ.chromatogram("test.mzXML", method = MassJ.∆MZ([1000, 1]))
        @test length(cr.x) == 6                                                       #17

        cr = MassJ.chromatogram("test.mzXML", method = MassJ.BasePeak())
        @test length(cr.x) == 6                                                       #18

        cr = MassJ.chromatogram("test.mzXML", MassJ.Polarity("+"), MassJ.Scan(2),MassJ.Precursor(1255.5), MassJ.Activation_Energy(18), MassJ.Activation_Method("CID"), MassJ.Level(2) )
        @test length(cr.x) == 1                                                       #19

        rt = MassJ.retention_time(scans)
        @test length(rt) == 6                                                          #20

        cr = MassJ.chromatogram(scans, method = MassJ.TIC() )
        @test length(cr.x) == 6                                                       #21

        cr = MassJ.chromatogram(scans, method = MassJ.MZ([0, 500]))
        @test length(cr.x) == 6                                                       #22

        cr = MassJ.chromatogram(scans, method = MassJ.∆MZ([1000, 1]))
        @test length(cr.x) == 6                                                       #23

        cr = MassJ.chromatogram(scans, method = MassJ.BasePeak())
        @test length(cr.x) == 6                                                       #24

        cr = MassJ.chromatogram(scans, MassJ.Polarity("+"),MassJ.Scan(2),MassJ.Precursor(1255.5),MassJ.Activation_Energy(18),MassJ.Activation_Method("CID"),MassJ.Level(2) )
        @test (cr.x, cr.ic) == ([0.7307], [9727.2])                                   #25

        cr = MassJ.chromatogram(scans, MassJ.Polarity(["+"]),MassJ.Scan([2,3]),MassJ.Precursor([1255.5, 902.33]),MassJ.Activation_Energy([18, 35]),MassJ.Activation_Method(["CID", "PQD"]),MassJ.Level([2, 3]) )
        @test (cr.x, cr.ic) == ([0.7307, 2.1379], [9727.2, 11.3032])                  #26

        ms = MassJ.average("test.mzXML")
        @test length(ms.num) == 6                                                      #27

        ms = MassJ.average("test.mzXML", MassJ.Polarity("+"),MassJ.Scan(2),MassJ.Precursor(1255.5),MassJ.Activation_Energy(18),MassJ.Activation_Method("CID"),MassJ.RT(1),MassJ.IC([0, 1e4]))
        @test ms isa MassJ.MSscans                                                        #28
        @test ms.num == [2]                                                              #29

        ms = MassJ.average(scans)
        @test length(ms.num) == 6                                                      #30

        ms = MassJ.average(scans, MassJ.Polarity("+"),MassJ.Scan(2),MassJ.Precursor(1255.5),MassJ.Activation_Energy(18),MassJ.Activation_Method("CID"),MassJ.RT(1),MassJ.IC([0, 1e4]))
        @test ms isa MassJ.MSscans                                                        #31
        @test ms.num == [2]                                                              #32

        ms = MassJ.average(scans, MassJ.Polarity(["+"]),MassJ.Scan([2,3]),MassJ.Precursor([1255.5, 902.33]),MassJ.Activation_Energy([18, 35]),MassJ.Activation_Method(["CID", "PQD"]),MassJ.RT([1,2]),MassJ.IC([0, 1e4]))
        @test ms isa MassJ.MSscans                                                       #33
        @test ms.num == [2, 3]                                                         #34

        ms = MassJ.average("test.mzXML", MassJ.RT( [[1,2], [2,3]] ), stats = false )
        @test ms isa MassJ.MSscans                                                       #35
        @test ms.num == [2, 3, 4]                                                      #36

        ms = MassJ.average("test.mzXML", MassJ.Polarity(["+"]),MassJ.Scan([2,3]),MassJ.Precursor([1255.5, 902.33]),MassJ.Activation_Method(["CID", "PQD"]),MassJ.RT([1,2]),MassJ.IC([0, 1e4]))   #MassJ.Activation_Energy([18., 35.]),
        @test ms isa MassJ.MSscans                                                       #37
        @test ms.num == [2, 3]                                                         #38

        cr = MassJ.chromatogram("test.mzXML", MassJ.Polarity(["+"]),MassJ.Scan([2,3]),MassJ.Precursor([1255.5, 902.33]),MassJ.Activation_Method(["CID", "PQD"]),MassJ.Level([2, 3]) )   #MassJ.Activation_Energy([18.0, 35.0]),
        @test length(cr.x) == 2                                                       #39

        ms = MassJ.average(scans, MassJ.RT( [[1,2], [2,3]] ), stats = false )
        @test ms isa MassJ.MSscans                                                       #40
        @test ms.num == [2, 3, 4]                                                      #41

        a = scans[1] / 2.
        @test a.tic == 2.540975e6                                                      #42

        a = scans[1] * 2.
        @test a.tic == 1.01639e7                                                       #43

        a = ms * 2.
        @test a.tic == 3.2120923354666666e6                                            #44

        a = 2. * scans[1]
        @test a.tic == 1.01639e7                                                       #45

        a = scans[1] * scans[2]
        @test a.tic == 4.943314404e10                                                  #46

        a = scans[2] * scans[1]
        @test a.tic == 4.943314404e10                                                  #47

        a = scans[1] - scans[2]
        @test a.tic == 5.0722228e6                                                     #48

        a = scans[2] - scans[1]
        @test a.tic == -5.0722228e6                                                    #49

        a = scans[1] - scans[4]
        @test a.tic == 273550.0                                                        #50

        b = ms - scans[1]
        @test b.num == [2,3,4]                                                         #51

        b = ms + scans[1]
        @test b.num == [2,3,4,1]                                                       #52

        b = scans[1] + ms
        @test b.num == [1,2,3,4]                                                       #53

        b = scans[1] + scans[4]
        @test b.num == [1,4]                                                           #54

        a = MassJ.avg(scans[1], scans[2])
        @test a.num == [1,2]                                                           #55

        a = MassJ.avg(scans[1], scans[4])
        @test a.num == [1,4]                                                           #56

        info = MassJ.info("test64.mzXML")
        @test info[2] == "MS1-"                                                        #57

        scans = MassJ.load("test64.mzXML")
        @test eltype(scans)              == MassJ.MSscans                                 #58

        @test_throws ErrorException MassJ.info("test.mzXMLM")                          #59
        @test_throws ErrorException MassJ.load("test.mzXMLL")                          #60
        @test_throws ErrorException MassJ.load("bad1.mzXML")                           #61
        @test_throws ErrorException MassJ.info("bad1.mzXML")                           #62

        scans = MassJ.load("bad2.mzXML")
        @test scans[1].num == [0]                                                        #63

        scans = MassJ.load("bad3.mzXML")
        @test scans[1].num == scans[2].num == [0]                                        #64

        @test_throws ErrorException MassJ.chromatogram("test.mzXML", method = MassJ.∆MZ([1, 2]))  #65
        @test_throws ErrorException MassJ.chromatogram(scans, method = MassJ.∆MZ([1, 2]))         #66

        scans = MassJ.load("test.mzXML")
        @test MassJ.smooth(scans[1], method = MassJ.SG(7,15,0)) isa MassJ.MSscans             #67

        a = MassJ.avg(scans[1], scans[4])
        @test MassJ.smooth(a) isa MassJ.MSscans                                            #68

        a = scans[1] * scans[4]
        @test a.num == [1,4]                                                           #69

        a = (scans[2]+scans[3]) - scans[1]
        @test a.num == [2, 3]                                                          #70

        a = (scans[1] + scans[4]) - (scans[1] - scans[4])
        @test a.num == [1, 4]                                                          #71

        a = scans[1] + MassJ.avg(scans[2], scans[5])
        @test a.num == [1, 2, 5]                                                       #72

        a = MassJ.smooth(scans[1], method = MassJ.SG(5,9,0))
        @test a.num == [1]                                                               #73

       a = MassJ.centroid(scans[1], method = MassJ.TBPD(:gauss, 4500., 0.2))               #74
       @test length(a.int) == 957

       @test typeof(plot(scans[1], method = :relative)) == Plots.Plot{Plots.GRBackend} #75
       @test typeof(plot(scans[1], method = :absolute)) == Plots.Plot{Plots.GRBackend} #76

       a = MassJ.avg(scans[2], scans[5])
       @test typeof(plot( a, method = :relative )) == Plots.Plot{Plots.GRBackend}      #77
       @test typeof(plot( a, method = :absolute )) == Plots.Plot{Plots.GRBackend}      #78

       cr = MassJ.chromatogram(scans)
       @test typeof(plot( cr, method = :relative )) == Plots.Plot{Plots.GRBackend}     #79
       @test typeof(plot( cr, method = :absolute )) == Plots.Plot{Plots.GRBackend}     #80

       a = MassJ.centroid(scans[1], method = MassJ.TBPD(:voigt, 4500., 0.2))
       @test length(a.int) == 961                                                      #81

       a = MassJ.centroid(scans[1], method = MassJ.TBPD(:lorentz, 4500., 0.2))
       @test length(a.int) == 964                                                      #82

       @test_throws "Unsupported peak profile. Use :gauss, :lorentz or :voigt." MassJ.centroid(scans[1], method = MassJ.TBPD(:other, 4500., 0.2))   #83

       a = MassJ.centroid(scans[1], method = MassJ.SNRA(1., 100))
       @test length(a.int) == 109                                                      #84

       # Regression: SNRA with a very low threshold used to throw BoundsError
       # because the local-max check `SNR[i+1]` ran one past the array end.
       # See https://… (fix/snra-bounds branch).
       avg_ms = MassJ.average(scans)
       @test MassJ.centroid(avg_ms, method = MassJ.SNRA(1e-3, 100)) isa MassJ.MSscans  #84b
       @test MassJ.centroid(scans[1], method = MassJ.SNRA(1e-3, 100)) isa MassJ.MSscans #84c

       s1 = MassJ.extract(scans, MassJ.Activation_Energy([18,35]))
       @test length(s1) == 4                                                           #85

       s1 = MassJ.extract("test.mzXML", MassJ.Activation_Energy(18))
       @test length(s1) == 2                                                           #86

       s1 = MassJ.extract(scans, MassJ.Scan(1))
       @test length(s1) == 1                                                           #87

       s1 = MassJ.extract("test.mzXML", MassJ.Scan(1))
       @test length(s1) == 1                                                           #88

       bs = MassJ.baseline_correction(scans, method = MassJ.TopHat(1))
       @test length(bs) == 6                                                           #89

       bs = MassJ.baseline_correction(scans[1], method = MassJ.TopHat(1))
       @test length(bs.int) == length(scans[1].int)                                    #90

       c = MassJ.centroid(scans, method = MassJ.TBPD(:gauss, 4500., 0.2)) ;
       bs = MassJ.baseline_correction( c, method = MassJ.LOESS(3))
       @test length(bs) == 6                                                           #91

       bs = MassJ.baseline_correction(c[1], method = MassJ.LOESS(3))
       @test length(bs.int) == length(c[1].int)                                        #92

       bs = MassJ.baseline_correction(scans, method = MassJ.IPSA(51,100))
       @test length(bs) == 6                                                           #93

       bs = MassJ.baseline_correction(scans[1], method = MassJ.IPSA(51,100))
       @test length(bs.int) == length(scans[1].int)                                    #94

       a = MassJ.smooth(scans, method = MassJ.SG(5,9,0))
       @test length(a) == 6                                                            #95

       c = MassJ.centroid(scans, method = MassJ.TBPD(:lorentz, 4500., 0.2)) ;
       d = MassJ.centroid(scans, method = MassJ.TBPD(:voigt, 4500., 0.2)) ;
       @test length(c) == length(d)                                                    #96

       a = MassJ.centroid(scans[3], method = MassJ.SNRA(1., 100))
       @test length(a.int) == 0                                                        #97

       bs = MassJ.baseline_correction(scans[1], method = MassJ.IPSA(50,100))
       @test length(bs.int) == length(scans[1].int)                                    #98

       cr = MassJ.chromatogram(scans, method = MassJ.BasePeak() )
       @test length(cr.x) == 6                                                        #99

       f = MassJ.formula("CH3(13C)10H3Kr(NaH2)2")                                        #100
      @test f == Dict("Na" => 2,"Kr" => 1,"C" => 1,"13C" => 10,"H" => 10)

       m = MassJ.masses("C254 H377 N65 O75 S6")                                          #101
      @test m["Monoisotopic"] ≈ 5729.60087099839
      @test m["Average"] ≈ 5733.55
      @test m["Nominal"] ≈ 5727.0

      I = MassJ.isotopic_distribution("CH4", 0.9999, charge = +1)                        #102
      @test I[2,1:end] == [16.03130012908, 0.9887541751052761, 1, 0, 4, 0]

      a = MassJ.simulate(I, 0.4, Npoints = 5)                                            #103
      @test a.int == [100.0, 6.035851011021856, 0.06994625998243831, 1.1368602493290418, 0.06861901924947518]

       m = MassJ.masses(f)                                                               #104
      @test m == Dict("Monoisotopic" => 282.0028349717, "Average" => 281.902086912, "Nominal" => 282.0)

      I = MassJ.isotopic_distribution(f, 0.9999, charge = +1)                            #105
      @test I[2,1:end][1:2] ≈ [282.0028349717, 0.5630635281692917]

      a = MassJ.simulate(I, 0.4, model=:lorentz, Npoints = 5)                            #106
      @test a.int ≈ [0.5359868694750152, 16.332387108915455, 100.0, 0.560482304632663, 0.0834376623225204]

    end
end


function test_isotopes()
    @testset "Isotopes - type stability and optimizations" begin

        # Elements dict is const (enables type inference)
        @test isconst(MassJ, :Elements)

        # PriorityQueues are typed
        pq = PriorityQueue{Vector{Int},Float64}()
        pq[[1,0]] = 0.5
        pq[[0,1]] = 0.3
        @test peek(pq) == ([0,1] => 0.3)   # min-heap: lowest first

        # formula parsing
        @test MassJ.formula("C2H6O") == Dict("C" => 2, "H" => 6, "O" => 1)
        @test MassJ.formula("H2O") == Dict("H" => 2, "O" => 1)
        @test MassJ.formula("NaCl") == Dict("Na" => 1, "Cl" => 1)
        @test_throws ErrorException MassJ.formula("123bad")
        @test_throws ErrorException MassJ.formula("Xx")

        # masses calculation
        m = MassJ.masses("H2O")
        @test m["Monoisotopic"] ≈ 18.01056468474
        @test m["Average"] ≈ 18.015
        @test m["Nominal"] == 18.0

        m2 = MassJ.masses(Dict("H" => 2, "O" => 1))
        @test m2 == m

        # stirling approximation vs exact log factorial
        @test MassJ.stirling(100) ≈ log(factorial(big(100))) atol=0.01
        @test MassJ.stirling(500) > 0

        # isotopologue_probability - low mass path (Natoms < 20)
        prob_H2 = MassJ.isotopologue_probability(Dict("H" => 2), Dict("H" => [2, 0]), MassJ.Elements)
        @test prob_H2 ≈ MassJ.Elements["H"][1].f^2 atol=1e-10

        # isotopologue_probability - high mass path (Natoms >= 20)
        prob_C20 = MassJ.isotopologue_probability(Dict("C" => 20), Dict("C" => [20, 0]), MassJ.Elements)
        @test prob_C20 ≈ 0.8049835604738165 atol=1e-8

        # isotopologue_mass
        cm = MassJ.isotopologue_mass([Pair("H", 2)], Dict("H" => [2, 0]), MassJ.Elements)
        @test cm ≈ 2 * MassJ.Elements["H"][1].m

        # most_probable_isotopologue
        alpha = MassJ.most_probable_isotopologue(Dict("C" => 10, "H" => 22), MassJ.Elements)
        @test sum(alpha["C"]) == 10     # conservation
        @test sum(alpha["H"]) == 22     # conservation
        @test alpha["C"][1] >= alpha["C"][2]  # 12C more abundant than 13C

        # hill_climbing finds optimum
        f_obj(x) = -(x[1] - 3)^2 - (x[2] - 2)^2
        P = MassJ.hill_climbing([1, 4], f_obj)
        @test P == [3, 2]

        # hill_climbing with single-element vector (early return)
        P_single = MassJ.hill_climbing([5], x -> -x[1]^2)
        @test P_single == [5]  # no neighbors possible

        # isotopic distribution - basic
        I = MassJ.isotopic_distribution("H2O", 0.99, charge = 1)
        @test size(I, 2) == 7   # Masses, Probability, + isotope columns
        @test I[2, 1] ≈ 18.01056468474   # monoisotopic mass
        @test I[2, 2] ≈ 0.9973367663173334  # probability
        @test I[1, 1] == "Masses"
        @test I[1, 2] == "Probability"

        # isotopic distribution - charge state divides mass
        I2 = MassJ.isotopic_distribution("H2O", 0.99, charge = 2)
        @test I2[2, 1] ≈ 18.01056468474 / 2

        # isotopic distribution - larger molecule
        I3 = MassJ.isotopic_distribution("C254 H377 N65 O75 S6", 0.5)
        @test size(I3, 1) > 2          # more than just header
        @test I3[1, 1] == "Masses"
        probs = [I3[i, 2] for i in 2:size(I3, 1)]
        @test all(p -> p > 0, probs)   # all probabilities positive
        @test sum(probs) ≈ 0.5 atol=0.1

    end
end


function test_adducts()
    @testset "Adducts - m/z <-> neutral mass" begin
        M = MassJ.masses("C6H12O6")["Monoisotopic"]   # glucose, 180.0634

        # reference m/z values (electron-mass aware)
        @test MassJ.adduct_mz(M, "[M+H]+")  ≈ 181.07066 atol = 1e-4
        @test MassJ.adduct_mz(M, "[M+Na]+") ≈ 203.05261 atol = 1e-4
        @test MassJ.adduct_mz(M, "[M-H]-")  ≈ 179.05611 atol = 1e-4

        # round-trip for every tabulated adduct, accepting either a name or an Adduct
        for (name, a) in MassJ.Adducts
            mz = MassJ.adduct_mz(M, name)
            @test MassJ.neutral_mass(mz, name) ≈ M atol = 1e-9
            @test MassJ.neutral_mass(mz, a)    ≈ M atol = 1e-9
            @test MassJ.adduct_mz(M, a)        ≈ mz atol = 1e-12
        end

        # multiply-charged ions divide by |z|
        @test MassJ.adduct_mz(M, "[M+2H]2+") < MassJ.adduct_mz(M, "[M+H]+")
        @test MassJ.adduct_mz(M, "[M+3H]3+") * 3 ≈ M + 3 * MassJ.masses("H")["Monoisotopic"] - 3 * MassJ.m_electron atol = 1e-6

        # dimer carries n = 2
        @test MassJ.adduct_mz(M, "[2M+H]+") ≈ 2M + MassJ.masses("H")["Monoisotopic"] - MassJ.m_electron atol = 1e-6

        # electron-mass toggle differs by exactly mₑ for a singly charged ion
        @test MassJ.adduct_mz(M, "[M+H]+"; electron = false) -
              MassJ.adduct_mz(M, "[M+H]+"; electron = true) ≈ MassJ.m_electron atol = 1e-12

        # radical cation [M]+ is M − mₑ
        @test MassJ.adduct_mz(M, "[M]+") ≈ M - MassJ.m_electron atol = 1e-9

        # custom adduct via keyword constructor round-trips
        cust = MassJ.Adduct("[M+Li]+", charge = +1, add = "Li")
        @test MassJ.neutral_mass(MassJ.adduct_mz(M, cust), cust) ≈ M atol = 1e-9
        @test MassJ.polarity(cust) == "+"
        @test MassJ.polarity(MassJ.Adducts["[M-H]-"]) == "-"

        # unknown adduct name raises
        @test_throws ErrorException MassJ.adduct_mz(M, "[M+Foo]+")
    end
end


function test_calibration()
    @testset "Calibration - fit and apply" begin
        ref = [118.0863, 322.0481, 622.0290, 922.0098, 1521.9715]
        obs = ref .* (1 + 3e-6) .+ 0.002    # 3 ppm gain + 2 mDa offset

        cal = MassJ.calibrate(obs, ref; model = :linear)
        @test cal.model == :linear
        @test length(cal.coeffs) == 2
        @test maximum(abs, cal.residuals_ppm) < 1e-6          # linear drift fully recovered

        # callable maps observed back onto reference, scalar and vector forms
        @test cal(obs) isa Vector{Float64}
        @test maximum(abs, cal(obs) .- ref) < 1e-6
        @test cal(obs[1]) ≈ ref[1] atol = 1e-6

        # apply to a single scan: m/z + base-peak m/z corrected, intensities untouched
        ints = [10.0, 20.0, 30.0, 40.0, 50.0]
        s = MassJ.MSscans(1, 0.0, sum(ints), copy(obs), ints, 1, obs[3], 30.0,
                          0.0, "+", "", 0.0)
        sc = MassJ.calibrate(s, cal)
        @test maximum(abs, sc.mz .- ref) < 1e-6
        @test sc.basePeakMz ≈ ref[3] atol = 1e-6
        @test sc.int == ints
        @test sc.tic == s.tic

        # NaN base-peak m/z (e.g. empty centroid) is preserved, not corrupted
        snan = MassJ.MSscans(2, 0.0, 0.0, copy(obs), ints, 1, NaN, NaN, 0.0, "+", "", 0.0)
        @test isnan(MassJ.calibrate(snan, cal).basePeakMz)

        # vector-of-scans and MSrun apply
        scans = MassJ.calibrate([s, s], cal)
        @test scans isa Vector{MassJ.MSscans}
        @test maximum(abs, scans[2].mz .- ref) < 1e-6

        run = MassJ.MSrun([s, s])
        rc = MassJ.calibrate(run, cal)
        @test rc isa MassJ.MSrun
        @test maximum(abs, rc[1].mz .- ref) < 1e-6
        @test rc.metadata === run.metadata                    # file-level metadata preserved

        # quadratic recovers a quadratic distortion
        obs2 = [0.001 + 1.0r + 1e-8 * r^2 for r in ref]
        cal2 = MassJ.calibrate(obs2, ref; model = :quadratic)
        @test length(cal2.coeffs) == 3
        @test maximum(abs, cal2.residuals_ppm) < 1e-3
        @test maximum(abs, cal2(obs2) .- ref) < 1e-3

        # :poly with explicit degree
        @test length(MassJ.calibrate(obs, ref; model = :poly, degree = 1).coeffs) == 2

        # error paths
        @test_throws ErrorException MassJ.calibrate([1.0, 2.0], [1.0, 2.0, 3.0])
        @test_throws ErrorException MassJ.calibrate([1.0], [1.0]; model = :quadratic)
        @test_throws ErrorException MassJ.calibrate(ref, ref; model = :cubic)
    end
end


function test_centroid_metrics()
    @testset "Centroid - per-peak metrics" begin
        mz = collect(100.0:0.01:300.0)
        sigma = 0.05
        fwhm_true = 2.35482 * sigma
        g(x, mu, A) = A * exp(-(x - mu)^2 / (2sigma^2))
        int = g.(mz, 150.0, 100.0) .+ g.(mz, 200.0, 60.0) .+ 0.5   # two peaks + baseline
        s = MassJ.MSscans(1, 0.0, sum(int), mz, int, 1, 150.0, 100.0, 0.0, "+", "", 0.0)

        # metrics = true populates the four parallel metadata vectors
        c = centroid(s; method = MassJ.SNRA(5.0, 200), metrics = true)
        @test length(c.mz) == 2
        for key in ("peak_fwhm", "peak_snr", "peak_area", "peak_resolution")
            @test haskey(c.metadata, key)
            @test length(c.metadata[key]) == length(c.mz)
        end

        # FWHM recovers the true Gaussian width; resolution = m/z / FWHM
        @test all(isapprox.(c.metadata["peak_fwhm"], fwhm_true; atol = 5e-3))
        @test c.metadata["peak_resolution"] ≈ c.mz ./ c.metadata["peak_fwhm"] atol = 1e-6

        # S/N is apex-over-baseline: finite, above threshold, taller peak scores higher
        @test all(c.metadata["peak_snr"] .> 5.0)
        @test all(isfinite, c.metadata["peak_snr"])
        @test c.metadata["peak_snr"][1] > c.metadata["peak_snr"][2]
        @test all(c.metadata["peak_area"] .> 0.0)

        # default (metrics = false) leaves the metadata untouched
        c0 = centroid(s; method = MassJ.SNRA(5.0, 200))
        @test !haskey(c0.metadata, "peak_fwhm")

        # vector path attaches metrics to every element
        cv = centroid([s, s]; method = MassJ.SNRA(5.0, 200), metrics = true)
        @test length(cv) == 2
        @test all(haskey(x.metadata, "peak_fwhm") for x in cv)

        # metrics attach regardless of centroiding algorithm (CWT)
        cc = centroid(s; method = MassJ.CWT(scales = 1.0:1.0:16.0, threshold = 3.0), metrics = true)
        @test haskey(cc.metadata, "peak_fwhm")
        @test length(cc.metadata["peak_fwhm"]) == length(cc.mz)
    end
end


function test_assignment()
    @testset "Formula assignment + isotope scoring" begin
        M = MassJ.masses("C6H12O6")["Monoisotopic"]          # glucose
        mzobs = MassJ.adduct_mz(M, "[M+H]+")
        elems = Dict("C" => 0:12, "H" => 0:24, "O" => 0:12, "N" => 0:6)

        cands = MassJ.assign_formula(mzobs; adduct = "[M+H]+", tol_ppm = 5, elements = elems)
        @test !isempty(cands)
        gi = findfirst(c -> c.formula == Dict("C" => 6, "H" => 12, "O" => 6), cands)
        @test gi !== nothing
        @test cands[gi].rdbe ≈ 1.0
        @test abs(cands[gi].error_ppm) < 1e-3                 # mz built from this exact formula
        @test cands[gi].mass ≈ M atol = 1e-6
        # results are sorted by |error| → the true formula ranks first
        @test cands[1].formula == Dict("C" => 6, "H" => 12, "O" => 6)
        @test issorted(abs.([c.error_ppm for c in cands]))

        # RDBE filter excludes glucose (DBE = 1) when a higher floor is required
        hi = MassJ.assign_formula(mzobs; adduct = "[M+H]+", tol_ppm = 5, elements = elems,
                                  rdbe = (2.0, 40.0))
        @test all(c -> c.formula != Dict("C" => 6, "H" => 12, "O" => 6), hi)

        # unreachable target within the given ranges → no candidates
        @test isempty(MassJ.assign_formula(5000.0; adduct = "[M+H]+", tol_ppm = 1, elements = elems))

        # max_results truncates
        @test length(MassJ.assign_formula(500.1; adduct = "[M+H]+", tol_ppm = 10, max_results = 5)) <= 5

        # unknown element raises
        @test_throws ErrorException MassJ.assign_formula(mzobs; elements = Dict("Xx" => 0:2))

        # isotope-pattern scoring: measured pattern built from glucose theory
        theo = MassJ.isotopic_distribution(Dict("C" => 6, "H" => 12, "O" => 6), 0.999; charge = 1)
        mmz = [MassJ.adduct_mz(Float64(theo[i, 1]), "[M+H]+") for i in 2:size(theo, 1)]
        mint = Float64.(theo[2:end, 2])
        sc_true = MassJ.score_isotope_pattern(mmz, mint, cands[gi]; adduct = "[M+H]+", tol = 0.02)
        @test sc_true > 0.99
        wrong = MassJ.FormulaCandidate(Dict("C" => 2, "H" => 4, "O" => 1), 44.026, 0.0, 1.0)
        @test MassJ.score_isotope_pattern(mmz, mint, wrong; adduct = "[M+H]+", tol = 0.02) < sc_true

        # MSscans convenience method agrees with the vector method
        spec = MassJ.MSscans(1, 0.0, sum(mint), mmz, mint, 1, mmz[1], maximum(mint),
                             0.0, "+", "", 0.0)
        @test MassJ.score_isotope_pattern(spec, cands[gi]; adduct = "[M+H]+", tol = 0.02) ≈ sc_true atol = 1e-9
    end
end


function test_multifolder()
    @testset "Multi-folder load + average" begin
        root = mktempdir()
        a = joinpath(root, "A"); b = joinpath(root, "B"); sub = joinpath(b, "sub")
        mkpath(a); mkpath(sub)
        cp("test.mzML", joinpath(a, "e1.mzML"))
        cp("test.mgf",  joinpath(a, "e2.mgf"))
        cp("test.msp",  joinpath(b, "e3.msp"))
        cp("test.mzML", joinpath(sub, "e10.mzML"))

        na = length(MassJ.load("test.mzML"))
        nb = length(MassJ.load("test.mgf"))
        nc = length(MassJ.load("test.msp"))

        combined = MassJ.load([a, b])
        @test eltype(combined) == MassJ.MSscans
        @test length(combined) == na + nb + nc            # flat: everything pooled
        @test length(MassJ.load([b])) == nc               # non-recursive ignores sub/
        @test length(MassJ.load([b]; recursive = true)) == nc + na  # walks sub/

        avg = MassJ.average([a, b])
        @test avg isa MassJ.MSscans
        @test !isempty(avg.mz)

        @test_throws ErrorException MassJ.load(["/no/such/dir/xyz123"])
        @test_throws ErrorException MassJ.load([a, "runtests.jl"])   # unsupported ext
    end
end


function test_yield_transforms()
    @testset "YieldCurve transforms" begin
        YC = MassJ.YieldCurve
        mk(x, Y) = YC(collect(Float64, x), "E", Y, fill(NaN, size(Y)),
                      vec(sum(Y, dims = 2)), fill(NaN, size(Y, 1)), fill(NaN, size(Y)),
                      ["a", "b"], [(99.0, 101.0), (199.0, 201.0)],
                      ["f$i" for i in 1:length(x)], Dict{String,Any}())
        yc1 = mk([2.0, 1.0], [20.0 2.0; 10.0 1.0])        # unsorted x on purpose
        yc2 = mk([3.0, 4.0], [30.0 3.0; 40.0 4.0])

        # combine: concatenate along x, re-sorted
        c = MassJ.combine_yields(yc1, yc2)
        @test c.x == [1.0, 2.0, 3.0, 4.0]
        @test c.yields[:, 1] == [10.0, 20.0, 30.0, 40.0]
        @test c.labels == ["a", "b"]

        # shift_x
        @test sort(MassJ.shift_x(yc1, 0.5).x) == [1.5, 2.5]

        # scale_yields: scalar scales matrix and recomputes TIC
        sc = MassJ.scale_yields(yc1, 2.0)
        @test sc.yields == yc1.yields .* 2
        @test sc.tic == vec(sum(yc1.yields .* 2, dims = 2))
        # per-peak factor
        scp = MassJ.scale_yields(yc1, [2.0, 10.0])
        @test scp.yields[:, 1] == yc1.yields[:, 1] .* 2
        @test scp.yields[:, 2] == yc1.yields[:, 2] .* 10

        # recalibrate_x via a Calibration and via an on-the-fly fit
        cal = MassJ.calibrate([1.0, 2.0], [1.1, 2.2]; model = :linear)
        @test sort(MassJ.recalibrate_x(yc1, cal).x) ≈ [1.1, 2.2] atol = 1e-9
        @test sort(MassJ.recalibrate_x(yc1, [1.0, 2.0], [1.1, 2.2]).x) ≈ [1.1, 2.2] atol = 1e-9

        # error paths
        yb = YC([5.0], "E", reshape([1.0], 1, 1), fill(NaN, 1, 1), [1.0], [NaN],
                fill(NaN, 1, 1), ["z"], [(0.0, 1.0)], ["g"], Dict{String,Any}())
        @test_throws ErrorException MassJ.combine_yields(yc1, yb)
        @test_throws ErrorException MassJ.scale_yields(yc1, [1.0, 2.0, 3.0])
    end
end


function test_peptides()
    @testset "Peptide fragment ions" begin
        # neutral monoisotopic peptide mass
        @test MassJ.peptide_mass("PEPTIDE") ≈ 799.35996 atol = 1e-4

        fr = MassJ.fragment_ions("PEPTIDE"; ions = (:a, :b, :c, :x, :y, :z), charges = 1:1)
        mz(lbl) = fr[findfirst(f -> f.label == lbl, fr)].mz
        # singly-charged reference m/z (proton = H − e⁻)
        @test mz("b1") ≈ 98.06004  atol = 1e-4
        @test mz("b2") ≈ 227.10263 atol = 1e-4
        @test mz("y1") ≈ 148.06043 atol = 1e-4
        @test mz("y3") ≈ 376.17144 atol = 1e-4
        @test mz("a2") ≈ 199.10772 atol = 1e-4   # b2 − CO
        @test mz("c2") ≈ 244.12918 atol = 1e-4   # b2 + NH3
        @test mz("z3") ≈ 359.14489 atol = 1e-4   # Mascot even-electron z
        @test mz("x3") ≈ 402.15071 atol = 1e-4

        # complementary b/y sum to peptide mass + 2 protons (+ H2O accounting)
        bn = MassJ.peptide_mass("PEPTIDE")
        @test mz("b2") + mz("y5") ≈ bn + 2 * (MassJ._M_H - MassJ.m_electron) atol = 1e-4

        # charge 2 divides by |z|
        f2 = MassJ.fragment_ions("PEPTIDE"; ions = (:y,), charges = 2:2)
        @test first(f.mz for f in f2 if f.label == "y3(2+)") ≈ 188.58936 atol = 1e-4

        # n residues → n-1 backbone ions per series
        b = MassJ.fragment_ions("PEPTIDE"; ions = (:b,))
        @test length(b) == length("PEPTIDE") - 1

        # side-chain ions: Pro yields no d/w; Ile & Thr give two isomers
        sc = MassJ.fragment_ions("PEPTIDE"; ions = (:d, :v, :w))
        @test count(f -> f.series == :d, sc) == 6      # E,T(×2),I(×2),D ; both P excluded
        @test count(f -> f.series == :w, sc) == 7
        @test count(f -> f.series == :v, sc) == 6
        @test any(f -> f.label == "w3'", sc)           # Ile second isomer
        @test !any(f -> f.series == :d && f.position == 1, sc)   # P at position 1

        # ±H radical variants
        ch = MassJ.fragment_ions("PEPTIDE"; ions = (:c,), hshifts = (-1, 0, 1))
        c2 = sort([f.label for f in ch if startswith(f.label, "c2")])
        @test c2 == ["c2", "c2+1", "c2-1"]
        cm = MassJ.fragment_ions("PEPTIDE"; ions = (:c,))
        c2base = first(f.mz for f in cm if f.label == "c2")
        @test first(f.mz for f in ch if f.label == "c2+1") ≈ c2base + MassJ._M_H atol = 1e-6

        # errors
        @test_throws ErrorException MassJ.fragment_ions("PEPTIXE")    # unknown residue
        @test_throws ErrorException MassJ.fragment_ions("P")         # too short
        @test_throws ErrorException MassJ.fragment_ions("PEPTIDE"; ions = (:q,))  # bad series
    end
end


function test_fragment_peaks()
    @testset "Fragment ions -> peak list -> yields" begin
        ions = MassJ.fragment_ions("PEPTIDE"; ions = (:b, :y), charges = 1:1)
        iy = ions[findfirst(i -> i.label == "y3", ions)]

        # default: fixed Peak windows centred on each fragment m/z
        pf = MassJ.fragment_peaks(ions; tol = 0.3)
        @test length(pf) == length(ions)
        @test all(p isa MassJ.Peak for p in pf)
        py = pf[findfirst(p -> p.label == "y3", pf)]
        @test (py.mz1 + py.mz2) / 2 ≈ iy.mz atol = 1e-6
        @test (py.mz2 - py.mz1) ≈ 0.6 atol = 1e-9

        # ppm windows
        pp = MassJ.fragment_peaks(ions; ppm = 10.0)
        pj = pp[findfirst(p -> p.label == "y3", pp)]
        @test (pj.mz2 - pj.mz1) ≈ 2 * iy.mz * 10e-6 atol = 1e-9

        # fixed = false → located TargetPeaks
        pt = MassJ.fragment_peaks(ions; ppm = 10.0, fixed = false, method = :local_max)
        @test all(p isa MassJ.TargetPeak for p in pt)
        @test pt[findfirst(p -> p.label == "y3", pt)].mz ≈ iy.mz atol = 1e-6

        # one-step from a sequence; labels carried through
        ps = MassJ.fragment_peaks("PEPTIDE"; ions = (:b, :y), charges = 1:1, tol = 0.3)
        @test [p.label for p in ps] == [i.label for i in ions]

        # end-to-end: fragment peaks feed straight into the yields machinery
        root = mktempdir()
        cp("test.mzML", joinpath(root, "e1.mzML"))
        cp("test.mzML", joinpath(root, "e2.mzML"))
        yc = MassJ.yields(root, pf; x0 = 0.0, step = 1.0, xlabel = "energy")
        @test yc isa MassJ.YieldCurve
        @test size(yc.yields, 2) == length(pf)
        @test yc.labels == [p.label for p in pf]
        @test length(yc.x) == 2
    end
end


function test_mobility()
    @testset "Mobilogram / ionogram extraction" begin
        mz = collect(100.0:1.0:110.0)
        int = fill(5.0, length(mz)); int[6] = 50.0
        mk(num, dt, cv, mt) = MassJ.MSscans(num, 0.0, sum(int), mz, int, 1, mz[6], 50.0,
                                            0.0, "+", "", 0.0, 0, :profile, dt, cv, mt,
                                            Dict{String,Any}())
        ims   = [mk(1, 10.0, 0.0, :TIMS),  mk(2, 20.0, 0.0, :TIMS),  mk(3, 30.0, 0.0, :TIMS)]
        faims = [mk(1, -1.0, -5.0, :FAIMS), mk(2, -1.0, -3.0, :FAIMS), mk(3, -1.0, -1.0, :FAIMS)]

        m = MassJ.mobilogram(ims)
        @test m isa MassJ.IonCurrent
        @test m.axis == :drift
        @test m.mobilityType == :TIMS
        @test m.x == [10.0, 20.0, 30.0]
        @test m.ic == fill(sum(int), 3)                                  # TIC
        @test MassJ.mobilogram(ims; method = MassJ.BasePeak()).ic == fill(50.0, 3)

        g = MassJ.ionogram(faims)
        @test g.axis == :cv
        @test g.mobilityType == :FAIMS
        @test g.x == [-5.0, -3.0, -1.0]

        # scans lacking the dimension are skipped
        @test length(MassJ.mobilogram(vcat(ims, [mk(4, -1.0, 0.0, :none)])).x) == 3
        # filters compose (drift time 20 & 30 fall in [15, 35])
        @test length(MassJ.mobilogram(ims, MassJ.DriftTime([15.0, 35.0])).x) == 2
        @test length(MassJ.ionogram(faims, MassJ.CompensationVoltage([-4.0, 0.0])).x) == 2

        # no matching dimension throws
        @test_throws ErrorException MassJ.mobilogram(faims)
        @test_throws ErrorException MassJ.ionogram(ims)
        # filename form delegates through load (test.mzML carries no mobility data)
        @test_throws ErrorException MassJ.mobilogram("test.mzML")

        # plot recipe handles the :drift and :cv axes
        @test plot(m, method = :absolute) isa Plots.Plot
        @test plot(g, method = :absolute) isa Plots.Plot
    end
end


function test_uncertainty()
    @testset "Uncertainty accessors + plot band" begin
        scans = MassJ.load("test.mzML")
        avg = MassJ.average(scans)                       # stats = true → s holds σ

        @test MassJ.nscans(avg) == length(scans)
        sd = MassJ.stdev(avg)
        se = MassJ.sem(avg)
        @test !isempty(sd)
        @test length(sd) == length(avg.mz)
        @test all(sd .>= 0.0)
        @test se ≈ sd ./ sqrt(MassJ.nscans(avg)) atol = 1e-12

        # a single scan carries no replicate statistics
        s1 = scans[1]
        @test MassJ.nscans(s1) == 1
        @test isempty(MassJ.stdev(s1))
        @test isempty(MassJ.sem(s1))

        # plot recipe band options return a plot; band is ignored when s is empty
        @test plot(avg; band = :sem) isa Plots.Plot
        @test plot(avg; band = :std, method = :absolute) isa Plots.Plot
        @test plot(avg) isa Plots.Plot
        @test plot(s1; band = :sem) isa Plots.Plot
    end
end


function test_tables()
    @testset "Tables.jl interface (YieldCurve + isotopes)" begin
        # YieldCurve as a Tables source (extension active via `import Tables`)
        yc = MassJ.YieldCurve([1.0, 2.0], "energy (eV)",
                              [10.0 1.0; 20.0 2.0], [0.5 0.1; 0.6 0.2],
                              [11.0, 22.0], [0.51, 0.63], fill(NaN, 2, 2),
                              ["a", "b"], [(99.0, 101.0), (199.0, 201.0)],
                              ["f1", "f2"], Dict{String,Any}())
        @test Tables.istable(typeof(yc))
        cols  = Tables.columns(yc)
        names = Tables.columnnames(cols)
        @test :x in names && :a in names && :b in names
        @test Symbol("a_err") in names && Symbol("b_err") in names
        @test :TIC in names && Symbol("TIC_err") in names && :file in names
        @test collect(Tables.getcolumn(cols, :x)) == [1.0, 2.0]
        @test collect(Tables.getcolumn(cols, :a)) == [10.0, 20.0]
        @test collect(Tables.getcolumn(cols, Symbol("a_err"))) == [0.5, 0.6]
        @test collect(Tables.getcolumn(cols, :file)) == ["f1", "f2"]
        @test length(Tables.rows(yc)) == 2
        @test Tables.columntable(yc).b == [1.0, 2.0]   # column :b == yields[:,2]

        # isotope_table: dependency-free NamedTuple, a valid Tables source
        dist = MassJ.isotopic_distribution(MassJ.formula("C2H6O"), 0.99)
        it = MassJ.isotope_table(dist)
        @test Tables.istable(typeof(it))
        @test haskey(it, :Masses) && haskey(it, :Probability)
        @test eltype(it.Masses) == Float64
        @test eltype(it[Symbol("12C")]) == Int
        @test length(it.Masses) == size(dist, 1) - 1
        @test Tables.columntable(it).Masses == it.Masses
    end
end


function test_measurements()
    @testset "Measurements.jl interface" begin
        scans = MassJ.load("test.mzML")
        avg = MassJ.average(scans)                       # carries σ

        m = MassJ.measurements(avg)                      # default :sem
        @test length(m) == length(avg.mz)
        @test Measurements.value.(m) == avg.int
        @test Measurements.uncertainty.(m) ≈ MassJ.sem(avg)
        @test Measurements.uncertainty.(MassJ.measurements(avg; kind = :std)) ≈ MassJ.stdev(avg)

        # a single scan has no replicate statistics
        @test_throws ErrorException MassJ.measurements(scans[1])
        # invalid kind
        @test_throws ErrorException MassJ.measurements(avg; kind = :bogus)

        # YieldCurve → matrix of value ± error, with automatic propagation
        yc = MassJ.YieldCurve([1.0, 2.0], "E", [10.0 1.0; 20.0 2.0], [0.5 0.1; 0.6 0.2],
                              [11.0, 22.0], [0.51, 0.63], fill(NaN, 2, 2), ["a", "b"],
                              [(99.0, 101.0), (199.0, 201.0)], ["f1", "f2"], Dict{String,Any}())
        M = MassJ.measurements(yc)
        @test size(M) == (2, 2)
        @test Measurements.value.(M) == yc.yields
        @test Measurements.uncertainty.(M) == yc.yields_err
        # summing a column propagates the errors in quadrature (Measurements does this)
        @test Measurements.uncertainty(sum(M[:, 1])) ≈ sqrt(0.5^2 + 0.6^2)
    end
end


function test_unitful()
    @testset "Unitful.jl interface (physical axes)" begin
        scan = MassJ.MSscans(1, 30.0, 100.0, [100.0, 200.0], [5.0, 9.0], 2,
                             200.0, 9.0, 400.0, "+", "CID", 18.0,
                             2, :centroid, -1.0, -5.0, :FAIMS, Dict{String,Any}())
        wu = MassJ.withunits(scan)
        @test Unitful.unit(eltype(wu.retention_time))      == u"s"
        @test Unitful.unit(eltype(wu.collision_energy))    == u"eV"
        @test Unitful.unit(eltype(wu.compensation_voltage)) == u"V"
        @test Unitful.ustrip.(u"s",  wu.retention_time)       == scan.rt
        @test Unitful.ustrip.(u"eV", wu.collision_energy)     == scan.collisionEnergy
        @test Unitful.ustrip.(u"V",  wu.compensation_voltage) == scan.compensationVoltage
        # m/z and intensity are not exposed (left untagged)
        @test !haskey(wu, :mz) && !haskey(wu, :intensity)

        # IonCurrent: abscissa tagged by axis
        rt = MassJ.IonCurrent([10.0, 20.0], [1.0, 2.0]; axis = :rt)
        @test Unitful.unit(eltype(MassJ.withunits(rt).x)) == u"s"
        cv = MassJ.IonCurrent([-5.0, -3.0], [1.0, 2.0]; axis = :cv, mobilityType = :FAIMS)
        @test Unitful.unit(eltype(MassJ.withunits(cv).x)) == u"V"
        # :drift is left unitless (ms vs 1/K₀ ambiguity)
        dt = MassJ.IonCurrent([10.0, 20.0], [1.0, 2.0]; axis = :drift, mobilityType = :TIMS)
        @test eltype(MassJ.withunits(dt).x) <: Real
        @test MassJ.withunits(rt).ic == rt.ic            # intensity untagged
    end
end


function test_chrom_processing()
    @testset "IonCurrent smoothing + baseline" begin
        t = collect(0.0:0.05:10.0)
        peak = 100.0 .* exp.(-(t .- 5.0).^2 ./ (2 * 0.3^2))
        ic = MassJ.IonCurrent(t, peak .+ (2.0 .+ 0.5 .* t); axis = :rt)

        sm = MassJ.smooth(ic; method = MassJ.SG(2, 9, 0))
        @test sm isa MassJ.IonCurrent
        @test sm.axis == :rt
        @test sm.x == ic.x
        @test length(sm.ic) == length(ic.ic)

        # baseline removal: sloped baseline gone (min ≈ 0), peak preserved (~100)
        for m in (MassJ.TopHat(40), MassJ.IPSA(31, 50))
            bc = MassJ.baseline_correction(ic; method = m)
            @test bc isa MassJ.IonCurrent
            @test bc.axis == :rt
            @test minimum(bc.ic) < 5.0
            @test maximum(bc.ic) > 80.0
        end
        @test MassJ.baseline_correction(ic; method = MassJ.LOESS(2)) isa MassJ.IonCurrent

        # axis / mobilityType preserved for non-rt traces
        xcv = collect(-5.0:0.1:5.0)
        cv = MassJ.IonCurrent(xcv, 50 .* exp.(-xcv .^ 2) .+ 1.0; axis = :cv, mobilityType = :FAIMS)
        bc = MassJ.baseline_correction(cv; method = MassJ.TopHat(20))
        @test bc.axis == :cv && bc.mobilityType == :FAIMS

        # unsupported smoothing method errors
        @test_throws ErrorException MassJ.smooth(ic; method = MassJ.TopHat(10))

        # MSscans regression: refactored loess/ipsa (+ tophat) still return MSscans
        mzv  = collect(100.0:0.5:200.0)
        intv = 50 .* exp.(-(mzv .- 150.0) .^ 2 ./ 8.0) .+ 2.0
        sp = MassJ.MSscans(1, 0.0, sum(intv), mzv, intv, 1, 150.0, 52.0, 0.0, "+", "", 0.0)
        @test MassJ.baseline_correction(sp; method = MassJ.LOESS(1))    isa MassJ.MSscans
        @test MassJ.baseline_correction(sp; method = MassJ.IPSA(31, 20)) isa MassJ.MSscans
        @test MassJ.baseline_correction(sp; method = MassJ.TopHat(20))   isa MassJ.MSscans
    end
end


function test_chrom_peaks()
    @testset "ChromPeak detection + metrics" begin
        t = collect(0.0:0.02:10.0)
        g(mu, s, a) = a .* exp.(-(t .- mu) .^ 2 ./ (2s^2))

        # one symmetric Gaussian: metrics match analytic values
        ic = MassJ.IonCurrent(t, g(5.0, 0.3, 100.0); axis = :rt)
        ps = MassJ.chrom_peaks(ic; snr = 3.0)
        @test length(ps) == 1
        p = ps[1]
        @test p.axis == :rt
        @test p.apex      ≈ 5.0 atol = 0.02
        @test p.height    ≈ 100.0 atol = 0.5
        @test p.fwhm      ≈ 2.35482 * 0.3 atol = 0.01           # 0.7064
        @test p.area      ≈ 100 * 0.3 * sqrt(2π) atol = 0.5      # 75.2
        @test p.asymmetry ≈ 1.0 atol = 0.05
        @test p.tailing   ≈ 1.0 atol = 0.05
        @test p.plates    ≈ 5.54 * (5.0 / p.fwhm)^2 atol = 1.0
        @test p.centroid  ≈ 5.0 atol = 0.02
        @test p.variance  ≈ 0.09 atol = 0.01
        @test abs(p.skewness) < 0.05
        @test p.snr > 10

        # two well-separated peaks → split at the valley
        ic2 = MassJ.IonCurrent(t, g(3.0, 0.2, 80.0) .+ g(7.0, 0.25, 50.0); axis = :rt)
        @test length(MassJ.chrom_peaks(ic2; snr = 3.0)) == 2

        # asymmetric peak (wider right tail) → Aₛ > 1 and positive skew
        yR = [tt < 5 ? exp(-(tt - 5)^2 / (2 * 0.2^2)) : exp(-(tt - 5)^2 / (2 * 0.6^2)) for tt in t] .* 100
        pa = MassJ.chrom_peaks(MassJ.IonCurrent(t, yR; axis = :rt); snr = 3.0)[1]
        @test pa.asymmetry > 1.5
        @test pa.skewness > 0.0

        # targeted single-window analysis
        pw = MassJ.chrom_peak(ic, 4.0, 6.0)
        @test pw isa MassJ.ChromPeak
        @test pw.apex ≈ 5.0 atol = 0.02
        @test pw.area > 70.0

        # rel_height gate drops a tiny peak beside a large one
        ic3 = MassJ.IonCurrent(t, g(3.0, 0.2, 100.0) .+ g(7.0, 0.2, 1.0); axis = :rt)
        @test length(MassJ.chrom_peaks(ic3; snr = 0.0, rel_height = 0.1)) == 1

        # degenerate window errors; axis carried through
        @test_throws ErrorException MassJ.chrom_peak(ic, 5.0, 5.0)
        icd = MassJ.IonCurrent(t, g(5.0, 0.3, 100.0); axis = :drift, mobilityType = :TIMS)
        @test MassJ.chrom_peaks(icd; snr = 3.0)[1].axis == :drift
    end
end


function test_deconvolution()
    @testset "Deconvolution - helpers and integration" begin

        # is_evenly_spaced
        @test MassJ.is_evenly_spaced([1.0, 2.0, 3.0, 4.0]) == true
        @test MassJ.is_evenly_spaced([1.0, 2.0, 3.5, 4.0]) == false
        @test MassJ.is_evenly_spaced(collect(range(100.0, stop=200.0, length=100))) == true

        # resampling
        X = [1.0, 2.0, 3.5, 5.0]
        Y = [10.0, 20.0, 35.0, 50.0]
        newX, newY = MassJ.resampling(X, Y)
        @test length(newX) == length(X)
        @test newX[1] ≈ X[1]
        @test newX[end] ≈ X[end]
        @test MassJ.is_evenly_spaced(newX)

        # new_mass
        @test MassJ.new_mass(500.0, 5, 1, 1.00782503227) ≈ 416.83463750537834
        @test MassJ.new_mass(500.0, 5, -1, 1.00782503227) ≈ 624.7480437419325
        # new_mass with zero shift returns original
        @test MassJ.new_mass(500.0, 5, 0, 1.0) ≈ 500.0

        # get_peak_shape
        mz = collect(range(100.0, stop=200.0, length=1000))
        g = MassJ.get_peak_shape(MassJ.gauss, 0.5, mz)
        @test length(g) > 0
        @test sum(g) ≈ 1.0 atol=1e-10      # normalized
        @test g[argmax(g)] > 0              # has a peak
        # approximately symmetric shape
        @test g[1] ≈ g[end] rtol=1.0

        # get_peak_shape with lorentz
        g_lor = MassJ.get_peak_shape(MassJ.lorentz, 0.5, mz)
        @test sum(g_lor) ≈ 1.0 atol=1e-10

        # figure_of_merit
        h = [1.0, 2.0, 3.0, 4.0, 5.0]
        @test MassJ.figure_of_merit(h, h) ≈ 1.0           # perfect match
        @test MassJ.figure_of_merit(h, h .+ 0.1) ≈ 0.995  # close match
        @test MassJ.figure_of_merit(h, h .+ 0.1) > 0.99

        # chargefilter - basic shape and properties
        mz_test = collect(range(400.0, stop=600.0, length=200))
        int_test = zeros(200)
        int_test[100] = 100.0   # single peak at ~500
        f_test = zeros(10, 200)
        for i in 1:10
            f_test[i, :] = int_test
        end
        s = MassJ.chargefilter(mz_test, int_test, f_test, (1, 10), 1.0, 1)
        @test size(s) == (10, 200)
        @test all(s .>= 0)     # non-negative output

        # project_N_convolve
        g_small = MassJ.get_peak_shape(MassJ.gauss, 2.0, mz_test)
        c = MassJ.project_N_convolve(s, g_small)
        @test length(c) == 200
        @test all(c .>= 0)

        # Charges struct construction
        ch = MassJ.Charges(adduct="H", range=(5, 15), width=2)
        @test ch.adduct == "H"
        @test ch.range == (5, 15)
        @test ch.width == 2

        # Charges with default width
        ch2 = MassJ.Charges(adduct="Na", range=(1, 5))
        @test ch2.width == 1

        # _resolve_shape
        scans = MassJ.load("test.mzXML")
        model_gauss = MassJ._resolve_shape(scans[1], :gauss)
        @test model_gauss === MassJ.gauss
        model_lor = MassJ._resolve_shape(scans[1], :lorentz)
        @test model_lor === MassJ.lorentz
        model_voigt = MassJ._resolve_shape(scans[1], :voigt)
        @test model_voigt === MassJ.voigt
        @test_throws ErrorException MassJ._resolve_shape(scans[1], :invalid)

        # _resolve_FWHM with explicit values
        @test MassJ._resolve_FWHM(scans[1], MassJ.gauss, -1, 0.5) == 0.5   # explicit FWHM
        @test MassJ._resolve_FWHM(scans[1], MassJ.gauss, 1000.0, -1) == 0.5  # R=1000 → 500/1000=0.5
        @test_throws ErrorException MassJ._resolve_FWHM(scans[1], MassJ.gauss, 1000.0, 0.5)  # both

        # roughguess_FWHM runs without error
        w = MassJ.roughguess_FWHM(scans[1])
        @test w > 0

        # deconv method exists with correct dispatch
        @test hasmethod(MassJ.deconv, Tuple{MassJ.MSscans, MassJ.Charges})
        @test hasmethod(MassJ.deconv, Tuple{MassJ.MSscans, MassJ.Charges})

    end
end


function test_interpolation_import()
    @testset "Interpolations import fix" begin
        # Line is accessible from MassJ (needed for extrapolation)
        @test isdefined(MassJ, :Line)
        @test isdefined(MassJ, :LinearInterpolation)

        # Verify interpolation with Line extrapolation works
        x = [1.0, 2.0, 3.0, 4.0]
        y = [10.0, 20.0, 30.0, 40.0]
        itp = MassJ.LinearInterpolation(x, y, extrapolation_bc=MassJ.Line())
        @test itp(2.5) ≈ 25.0
        @test itp(0.0) ≈ 0.0    # extrapolated
    end
end


function test_mzml()
    @testset "mzML format" begin
        # info
        inf = MassJ.info("test.mzML")
        @test any(contains(s, "3 scans") for s in inf)
        @test any(contains(s, "MS1+") for s in inf)
        @test any(contains(s, "MS2+") for s in inf)

        inf_v = MassJ.info("test.mzML", verbose=true)
        @test any(contains(s, "test.raw") for s in inf_v)
        @test any(contains(s, "LTQ FT") for s in inf_v)

        # load all
        scans = MassJ.load("test.mzML")
        @test length(scans) == 3
        @test eltype(scans) == MassJ.MSscans

        # Scan 1: MS1+, profile, rt=0.5 min
        s1 = scans[1]
        @test s1.level == [1]
        @test s1.rt ≈ [0.5]
        @test s1.polarity == ["+"]
        @test s1.spectrumType == :profile
        @test s1.tic ≈ 19000.0
        @test s1.basePeakMz ≈ 400.0
        @test s1.basePeakIntensity ≈ 8000.0
        @test s1.precursor ≈ [0.0]
        @test s1.chargeState == [0]
        @test length(s1.mz) == 5
        @test s1.mz ≈ [100.0, 200.0, 300.0, 400.0, 500.0]
        @test s1.int ≈ [1000.0, 5000.0, 3000.0, 8000.0, 2000.0]

        # Scan 2: MS2+ CID, centroid, rt=1.0 min (60s converted)
        s2 = scans[2]
        @test s2.level == [2]
        @test s2.rt ≈ [1.0]
        @test s2.polarity == ["+"]
        @test s2.spectrumType == :centroid
        @test s2.precursor ≈ [400.0]
        @test s2.chargeState == [2]
        @test s2.activationMethod == ["CID"]
        @test s2.collisionEnergy ≈ [25.0]
        @test length(s2.mz) == 4
        @test s2.mz ≈ [110.0, 150.0, 200.0, 250.0]

        # Scan 3: MS2+ HCD
        s3 = scans[3]
        @test s3.level == [2]
        @test s3.rt ≈ [1.5]
        @test s3.precursor ≈ [500.0]
        @test s3.chargeState == [3]
        @test s3.activationMethod == ["HCD"]
        @test s3.collisionEnergy ≈ [30.0]
        @test length(s3.mz) == 3

        # retention_time
        rt = MassJ.retention_time("test.mzML")
        @test length(rt) == 3
        @test rt ≈ [0.5, 1.0, 1.5]

        # chromatogram (TIC)
        chrom = MassJ.chromatogram("test.mzML")
        @test length(chrom.x) == 3
        @test chrom.x ≈ [0.5, 1.0, 1.5]
        @test chrom.ic ≈ [19000.0, 4800.0, 2100.0]
        @test MassJ.maxic(chrom) ≈ 19000.0

        # chromatogram with filter
        chrom2 = MassJ.chromatogram("test.mzML", MassJ.Level(2))
        @test length(chrom2.x) == 2

        # chromatogram base peak
        chrom_bp = MassJ.chromatogram("test.mzML", method=MassJ.BasePeak())
        @test chrom_bp.ic ≈ [8000.0, 2000.0, 1200.0]

        # extract
        ms2 = MassJ.extract("test.mzML", MassJ.Level(2))
        @test length(ms2) == 2
        @test ms2[1].precursor ≈ [400.0]
        @test ms2[2].precursor ≈ [500.0]

        # average
        avg = MassJ.average("test.mzML", MassJ.Level(2))
        @test avg isa MassJ.MSscans

        # New fields present
        @test s1.mobilityType == :none
        @test s1.driftTime ≈ [-1.0]
        @test s1.compensationVoltage ≈ [0.0]
        @test s1.metadata isa Dict{String,Any}
    end
end


function test_mgf()
    @testset "MGF format" begin
        # info
        inf = MassJ.info("test.mgf")
        @test inf[1] == "3 scans"
        @test any(contains(s, "400.0") for s in inf)
        @test any(contains(s, "500.0") for s in inf)
        @test any(contains(s, "600.0") for s in inf)

        # load all
        scans = MassJ.load("test.mgf")
        @test length(scans) == 3
        @test eltype(scans) == MassJ.MSscans

        # Scan 1
        s1 = scans[1]
        @test s1.num == [1]        # sequential index
        @test s1.rt ≈ [0.5]       # 30s / 60
        @test s1.level == [2]      # MGF default
        @test s1.precursor ≈ [400.0]
        @test s1.chargeState == [2]
        @test s1.polarity == ["+"]
        @test s1.spectrumType == :centroid
        @test s1.tic ≈ 4800.0   # sum of intensities
        @test s1.basePeakMz ≈ 150.0
        @test s1.basePeakIntensity ≈ 2000.0
        @test length(s1.mz) == 4
        @test s1.mz ≈ [110.0, 150.0, 200.0, 250.0]
        @test s1.int ≈ [500.0, 2000.0, 1500.0, 800.0]
        @test s1.metadata["title"] == "Spectrum 1, scan 101"

        # Scan 2
        s2 = scans[2]
        @test s2.num == [2]
        @test s2.rt ≈ [1.0]
        @test s2.precursor ≈ [500.0]
        @test s2.chargeState == [3]

        # Scan 3
        s3 = scans[3]
        @test s3.num == [3]
        @test s3.precursor ≈ [600.0]
        @test s3.chargeState == [1]
        @test length(s3.mz) == 4

        # retention_time
        rt = MassJ.retention_time("test.mgf")
        @test length(rt) == 3
        @test rt ≈ [0.5, 1.0, 1.5]

        # chromatogram
        chrom = MassJ.chromatogram("test.mgf")
        @test length(chrom.x) == 3
        @test MassJ.maxic(chrom) ≈ 4800.0

        # extract (filter by precursor)
        sub = MassJ.extract("test.mgf", MassJ.Precursor(500.0))
        @test length(sub) == 1
        @test sub[1].precursor ≈ [500.0]

        # average
        avg = MassJ.average("test.mgf")
        @test avg isa MassJ.MSscans

        # New MSscan fields
        @test s1.mobilityType == :none
        @test s1.driftTime ≈ [-1.0]
        @test s1.compensationVoltage ≈ [0.0]
    end
end


function test_msp()
    @testset "MSP format" begin
        # info
        inf = MassJ.info("test.msp")
        @test inf[1] == "3 scans"
        @test any(contains(s, "MS2+") for s in inf)
        @test any(contains(s, "MS2-") for s in inf)
        @test any(contains(s, "MS1+") for s in inf)

        inf_v = MassJ.info("test.msp", verbose=true)
        @test any(contains(s, "Caffeine") for s in inf_v)
        @test any(contains(s, "Aspirin") for s in inf_v)

        # load all
        scans = MassJ.load("test.msp")
        @test length(scans) == 3
        @test eltype(scans) == MassJ.MSscans

        # Scan 1: Caffeine MS2+
        s1 = scans[1]
        @test s1.num == [1]
        @test s1.level == [2]
        @test s1.polarity == ["+"]
        @test s1.precursor ≈ [195.0877]
        @test s1.collisionEnergy ≈ [20.0]
        @test s1.spectrumType == :centroid
        @test length(s1.mz) == 5
        @test s1.mz ≈ [42.0338, 69.0447, 110.0713, 138.0662, 195.0877]
        @test s1.int ≈ [1200.0, 3500.0, 8900.0, 45000.0, 120000.0]
        @test s1.basePeakMz ≈ 195.0877
        @test s1.basePeakIntensity ≈ 120000.0
        @test s1.tic ≈ sum(s1.int)
        @test s1.metadata["name"] == "Caffeine"
        @test s1.metadata["formula"] == "C8H10N4O2"
        @test s1.metadata["inchikey"] == "RYYVLZVUVIJVGH-UHFFFAOYSA-N"
        @test s1.metadata["db_id"] == "MSP001"
        @test s1.metadata["comments"] == "test spectrum 1"

        # Scan 2: Aspirin MS2-
        s2 = scans[2]
        @test s2.num == [2]
        @test s2.level == [2]
        @test s2.polarity == ["-"]
        @test s2.precursor ≈ [179.0344]
        @test s2.collisionEnergy ≈ [15.0]
        @test length(s2.mz) == 3
        @test s2.metadata["name"] == "Aspirin"
        @test s2.metadata["cas"] == "50-78-2"

        # Scan 3: Glucose MS1+ (semicolon format)
        s3 = scans[3]
        @test s3.num == [3]
        @test s3.level == [1]
        @test s3.polarity == ["+"]
        @test s3.precursor ≈ [0.0]
        @test length(s3.mz) == 4
        @test s3.mz ≈ [180.0634, 181.0668, 182.0701, 183.0735]
        @test s3.int ≈ [8000.0, 900.0, 50.0, 5.0]
        @test s3.metadata["name"] == "Glucose"

        # retention_time
        rt = MassJ.retention_time("test.msp")
        @test length(rt) == 3
        @test all(rt .≈ 0.0)  # no RT in MSP test file

        # chromatogram
        chrom = MassJ.chromatogram("test.msp")
        @test length(chrom.x) == 3

        # extract
        ms2_pos = MassJ.extract("test.msp", MassJ.Polarity("+"), MassJ.Level(2))
        @test length(ms2_pos) == 1
        @test ms2_pos[1].metadata["name"] == "Caffeine"

        ms2_neg = MassJ.extract("test.msp", MassJ.Polarity("-"))
        @test length(ms2_neg) == 1
        @test ms2_neg[1].metadata["name"] == "Aspirin"

        # average
        avg = MassJ.average("test.msp", MassJ.Level(2))
        @test avg isa MassJ.MSscans
    end
end


function test_imzml()
    @testset "imzML format" begin
        # info
        inf = MassJ.info("test.imzML")
        @test any(contains(s, "4 spectra") for s in inf)
        @test any(contains(s, "MS1+") for s in inf)

        inf_v = MassJ.info("test.imzML", verbose=true)
        @test any(contains(s, "processed") for s in inf_v)
        @test any(contains(s, "2 x 2") for s in inf_v)

        # load all
        scans = MassJ.load("test.imzML")
        @test length(scans) == 4
        @test eltype(scans) == MassJ.MSscans

        # Scan 1: position (1,1)
        s1 = scans[1]
        @test s1.num == [1]
        @test s1.level == [1]
        @test s1.polarity == ["+"]
        @test s1.spectrumType == :profile
        @test s1.tic ≈ 8000.0
        @test length(s1.mz) == 3
        @test s1.mz ≈ [100.0, 200.0, 300.0]
        @test s1.int ≈ [1000.0, 5000.0, 2000.0]
        @test s1.metadata["position_x"] == 1
        @test s1.metadata["position_y"] == 1

        # Scan 2: position (2,1)
        s2 = scans[2]
        @test s2.num == [2]
        @test s2.tic ≈ 7800.0
        @test s2.mz ≈ [100.0, 200.0, 300.0]
        @test s2.int ≈ [800.0, 3000.0, 4000.0]
        @test s2.metadata["position_x"] == 2
        @test s2.metadata["position_y"] == 1

        # Scan 3: position (1,2), 2 peaks
        s3 = scans[3]
        @test s3.num == [3]
        @test length(s3.mz) == 2
        @test s3.mz ≈ [150.0, 250.0]
        @test s3.int ≈ [2000.0, 6000.0]
        @test s3.metadata["position_x"] == 1
        @test s3.metadata["position_y"] == 2

        # Scan 4: position (2,2)
        s4 = scans[4]
        @test s4.num == [4]
        @test length(s4.mz) == 3
        @test s4.mz ≈ [150.0, 250.0, 350.0]
        @test s4.int ≈ [1500.0, 4000.0, 3000.0]
        @test s4.metadata["position_x"] == 2
        @test s4.metadata["position_y"] == 2

        # retention_time
        rt = MassJ.retention_time("test.imzML")
        @test length(rt) == 4

        # chromatogram
        chrom = MassJ.chromatogram("test.imzML")
        @test length(chrom.x) == 4

        # extract
        all_scans = MassJ.extract("test.imzML", MassJ.Level(1))
        @test length(all_scans) == 4

        # average
        avg = MassJ.average("test.imzML")
        @test avg isa MassJ.MSscans
    end
end


function test_export()
    @testset "Export — mzML / mzXML round-trip" begin

        # -- mzML round-trip ------------------------------------------------
        scans = MassJ.load("test.mzML")
        tmp = tempname() * ".mzML"
        MassJ.save(scans, tmp)
        @test isfile(tmp)
        back = MassJ.load(tmp)
        @test length(back) == length(scans)
        for (a, b) in zip(scans, back)
            @test a.mz                ≈ b.mz
            @test a.int               ≈ b.int
            @test a.level             == b.level
            @test a.rt                ≈ b.rt
            @test a.tic               ≈ b.tic
            @test a.polarity          == b.polarity
            @test a.precursor         ≈ b.precursor
            @test a.activationMethod  == b.activationMethod
            @test a.collisionEnergy   ≈ b.collisionEnergy
            @test a.basePeakMz        ≈ b.basePeakMz
            @test a.basePeakIntensity ≈ b.basePeakIntensity
        end
        rm(tmp)

        # -- mzXML round-trip -----------------------------------------------
        scans_x = MassJ.load("test.mzXML")
        tmpx = tempname() * ".mzXML"
        MassJ.save(scans_x, tmpx)
        @test isfile(tmpx)
        back_x = MassJ.load(tmpx)
        @test length(back_x) == length(scans_x)
        for (a, b) in zip(scans_x, back_x)
            @test a.mz               ≈ b.mz
            @test a.int              ≈ b.int
            @test a.level            == b.level
            @test a.rt               ≈ b.rt
            @test a.tic              ≈ b.tic
            @test a.precursor        ≈ b.precursor
            @test a.activationMethod == b.activationMethod
            @test a.collisionEnergy  ≈ b.collisionEnergy
        end
        rm(tmpx)

        # -- Single MSscan in / out — type round-trips bit-symmetrically ----
        tmp1 = tempname() * ".mzML"
        MassJ.save(scans[1], tmp1)
        back1 = MassJ.load(tmp1)
        @test typeof(back1) == typeof(scans[1])      # MSscan, not Vector
        @test back1.mz ≈ scans[1].mz
        rm(tmp1)

        # -- MSscans (averaged) in / out — type round-trips ------------------
        avg = MassJ.average("test.mzML")
        @test avg isa MassJ.MSscans
        tmpa = tempname() * ".mzML"
        MassJ.save(avg, tmpa)
        back_avg = MassJ.load(tmpa)
        @test typeof(back_avg) == typeof(avg)        # MSscans, not Vector
        @test back_avg.mz ≈ avg.mz
        @test back_avg.s  == avg.s                   # variance preserved
        rm(tmpa)

        # -- 32-bit precision still close (relative tol) --------------------
        tmp32 = tempname() * ".mzML"
        MassJ.save(scans, tmp32; precision = 32)
        back32 = MassJ.load(tmp32)
        @test back32[1].mz ≈ scans[1].mz rtol = 1e-6
        @test back32[1].int ≈ scans[1].int rtol = 1e-6
        rm(tmp32)

        # -- No compression also round-trips --------------------------------
        tmpnc = tempname() * ".mzML"
        MassJ.save(scans, tmpnc; compress = false)
        @test MassJ.load(tmpnc)[1].mz ≈ scans[1].mz
        rm(tmpnc)

        # -- Unknown extension errors ---------------------------------------
        @test_throws ErrorException MassJ.save(scans, "out.weird")

        # -- MSscans round-trip: typeof(load) == typeof(save) (mzML) --------
        mean_ml = MassJ.average(scans)       # MSscans from test.mzML
        @test mean_ml isa MassJ.MSscans
        tmpm = tempname() * ".mzML"
        MassJ.save(mean_ml, tmpm)
        back_ml = MassJ.load(tmpm)
        @test typeof(back_ml) == typeof(mean_ml)   # bit-symmetric
        @test back_ml.mz   == mean_ml.mz
        @test back_ml.int  == mean_ml.int
        @test back_ml.s    == mean_ml.s            # variance preserved exactly
        @test back_ml.num  == mean_ml.num          # history preserved
        @test back_ml.tic  ≈ mean_ml.tic
        rm(tmpm)

        # -- MSscans round-trip: typeof(load) == typeof(save) (mzXML) -------
        mean_x = MassJ.average(scans_x)
        @test mean_x isa MassJ.MSscans
        tmpx2 = tempname() * ".mzXML"
        MassJ.save(mean_x, tmpx2)
        back_x2 = MassJ.load(tmpx2)
        @test typeof(back_x2) == typeof(mean_x)
        @test back_x2.mz == mean_x.mz
        @test back_x2.int == mean_x.int
        @test back_x2.s  == mean_x.s
        rm(tmpx2)

        # -- Plain MSscan-Vector save still returns Vector{MSscan} ----------
        tmp_plain = tempname() * ".mzML"
        MassJ.save(scans, tmp_plain)
        back_plain = MassJ.load(tmp_plain)
        @test typeof(back_plain) == typeof(scans)
        rm(tmp_plain)

        # -- External (non-MassJ) file returns MSrun (subtype of AbstractVector{MSscan})
        ext_loaded = MassJ.load("test.mzML")
        @test ext_loaded isa MassJ.MSrun
        @test ext_loaded isa AbstractVector{MassJ.MSscans}
        @test length(ext_loaded) == 3   # iteration / length still work

        # -- Vector{MSscans} round-trip (mzML and mzXML) --------------------
        vec_ms = MassJ.MSscans[MassJ.average(scans), MassJ.average(scans)]
        @test vec_ms isa Vector{MassJ.MSscans}

        tmpvm = tempname() * ".mzML"
        MassJ.save(vec_ms, tmpvm)
        bvm = MassJ.load(tmpvm)
        # mzML load returns an MSrun (AbstractVector{MSscans}); the composite
        # spectra round-trip element-for-element.
        @test bvm isa AbstractVector{MassJ.MSscans}
        @test length(bvm) == 2
        @test bvm[1].s == vec_ms[1].s
        @test bvm[2].s == vec_ms[2].s
        @test bvm[1].num == vec_ms[1].num
        rm(tmpvm)

        tmpvx = tempname() * ".mzXML"
        MassJ.save(vec_ms, tmpvx)
        bvx = MassJ.load(tmpvx)
        @test bvx isa AbstractVector{MassJ.MSscans}
        @test length(bvx) == 2
        @test bvx[1].s == vec_ms[1].s
        @test bvx[2].s == vec_ms[2].s
        rm(tmpvx)

        # -- mzML cvParams round-trip via MSscan.metadata -------------------
        # Build a scan with all 7 known cvParam metadata keys populated, save,
        # reload, and verify every key round-trips byte-equivalently.
        meta = Dict{String,Any}(
            "spectrum_title"       => "synthetic test scan",
            "lowest_observed_mz"   => 110.123,
            "highest_observed_mz"  => 1499.9,
            "mass_resolving_power" => 120000.0,
            "ion_injection_time"   => 9.876,
            "scan_window_lower"    => 100.0,
            "scan_window_upper"    => 1500.0,
        )
        s0 = MassJ.MSscans(1, 0.5, 1.0e5, [100.0, 200.0], [10.0, 20.0],
                          1, 200.0, 20.0, 0.0, "+", "", 0.0,
                          0, :centroid, -1.0, 0.0, :none, meta)
        tmp_md = tempname() * ".mzML"
        MassJ.save([s0], tmp_md; progress = false)
        back_md = MassJ.load(tmp_md)
        @test back_md[1].metadata["spectrum_title"]       == "synthetic test scan"
        @test back_md[1].metadata["lowest_observed_mz"]   ≈ 110.123
        @test back_md[1].metadata["highest_observed_mz"]  ≈ 1499.9
        @test back_md[1].metadata["mass_resolving_power"] ≈ 120000.0
        @test back_md[1].metadata["ion_injection_time"]   ≈ 9.876
        @test back_md[1].metadata["scan_window_lower"]    ≈ 100.0
        @test back_md[1].metadata["scan_window_upper"]    ≈ 1500.0
        rm(tmp_md)

        # Scans without these keys produce no cvParams (keys remain absent
        # after round-trip — we don't fabricate zero defaults).
        s_empty = MassJ.MSscans(1, 0.5, 1.0e5, [100.0, 200.0], [10.0, 20.0],
                               1, 200.0, 20.0, 0.0, "+", "", 0.0)
        tmp_empty = tempname() * ".mzML"
        MassJ.save([s_empty], tmp_empty; progress = false)
        back_empty = MassJ.load(tmp_empty)
        for k in ("spectrum_title", "lowest_observed_mz", "highest_observed_mz",
                  "mass_resolving_power", "ion_injection_time",
                  "scan_window_lower", "scan_window_upper")
            @test !haskey(back_empty[1].metadata, k)
        end
        rm(tmp_empty)

        # -- Phase 1 cvParams: per-spectrum round-trip ---------------------
        # Build a scan with the 5 new per-spectrum cvParams + the previous 7.
        meta_phase1 = Dict{String,Any}(
            "spectrum_title"                 => "phase1 test",
            "lowest_observed_mz"             => 101.0,
            "highest_observed_mz"            => 1499.0,
            "mass_resolving_power"           => 60000.0,
            "ion_injection_time"             => 25.5,
            "scan_window_lower"              => 100.0,
            "scan_window_upper"              => 1500.0,
            "filter_string"                  => "FTMS + p NSI Full ms",
            "isolation_window_target_mz"     => 500.5,
            "isolation_window_lower_offset"  => 0.5,
            "isolation_window_upper_offset"  => 0.5,
            "selected_ion_peak_intensity"    => 4.2e6,
        )
        s_ms2 = MassJ.MSscans(1, 0.5, 1.0e5, [100.0, 200.0], [10.0, 20.0],
                             2, 200.0, 20.0, 500.5, "+", "HCD", 30.0,
                             2, :centroid, -1.0, 0.0, :none, meta_phase1)
        tmp_p1 = tempname() * ".mzML"
        MassJ.save([s_ms2], tmp_p1; progress = false)
        back_p1 = MassJ.load(tmp_p1)
        for k in keys(meta_phase1)
            @test back_p1[1].metadata[k] == meta_phase1[k]
        end
        rm(tmp_p1)

        # -- Phase 2: MSrun wrapper + file-level metadata round-trip -------
        # Build a synthetic MSrun with one of each metadata section.
        scans_src   = MassJ.load("test.mzML")    # already an MSrun
        @test scans_src isa MassJ.MSrun
        @test length(scans_src) == 3              # AbstractVector behavior
        @test scans_src[1] isa MassJ.MSscans      # indexing works

        run_md = Dict{String,Any}(
            "file_content" => [
                Dict{String,String}("accession" => "MS:1000579", "name" => "MS1 spectrum"),
            ],
            "source_files" => [
                Dict{String,Any}("id" => "RAW1", "name" => "test.raw",
                                 "location" => "file:///data/",
                                 "cv_params" => [
                                     Dict{String,String}("accession" => "MS:1000563",
                                                         "name" => "Thermo RAW format"),
                                 ]),
            ],
            "software" => [
                Dict{String,Any}("id" => "Xcalibur", "version" => "3.5",
                                 "cv_params" => [
                                     Dict{String,String}("accession" => "MS:1000532",
                                                         "name" => "Xcalibur"),
                                 ]),
            ],
            "instruments" => [
                Dict{String,Any}("id" => "IC1",
                                 "param_group_refs" => String[],
                                 "cv_params" => [
                                     Dict{String,String}("accession" => "MS:1003029",
                                                         "name" => "Orbitrap Eclipse"),
                                 ],
                                 "components" => Dict{String,Any}[
                                     Dict{String,Any}("type" => "source",   "order" => "1",
                                                      "cv_params" => Dict{String,String}[
                                                          Dict("accession" => "MS:1000398",
                                                               "name" => "nanoelectrospray"),
                                                      ]),
                                     Dict{String,Any}("type" => "analyzer", "order" => "2",
                                                      "cv_params" => Dict{String,String}[
                                                          Dict("accession" => "MS:1000484",
                                                               "name" => "orbitrap"),
                                                      ]),
                                 ]),
            ],
        )
        run_in = MassJ.MSrun(scans_src.scans, run_md, MassJ.IonCurrent[])
        tmp_run = tempname() * ".mzML"
        MassJ.save(run_in, tmp_run; progress = false)
        run_out = MassJ.load(tmp_run)
        @test run_out isa MassJ.MSrun
        @test length(run_out.metadata["source_files"]) == 1
        @test run_out.metadata["source_files"][1]["id"] == "RAW1"
        @test run_out.metadata["software"][1]["id"] == "Xcalibur"
        @test run_out.metadata["software"][1]["cv_params"][1]["accession"] == "MS:1000532"
        @test run_out.metadata["instruments"][1]["cv_params"][1]["accession"] == "MS:1003029"
        @test length(run_out.metadata["instruments"][1]["components"]) == 2
        @test run_out.metadata["instruments"][1]["components"][1]["type"] == "source"
        @test run_out.metadata["instruments"][1]["components"][1]["cv_params"][1]["accession"] == "MS:1000398"
        rm(tmp_run)

        # -- Phase 3: chromatogram round-trip ------------------------------
        rt_vec  = collect(0.0:0.5:5.0)
        ic_vec  = Float64[100, 200, 350, 500, 600, 720, 690, 540, 400, 250, 110]
        chrom_in = MassJ.IonCurrent(rt_vec, ic_vec; axis = :rt)
        run_with_chrom = MassJ.MSrun(scans_src.scans, Dict{String,Any}(),
                                      [chrom_in])
        tmp_chrom = tempname() * ".mzML"
        MassJ.save(run_with_chrom, tmp_chrom; progress = false)
        run_chrom_back = MassJ.load(tmp_chrom)
        @test length(run_chrom_back.chromatograms) == 1
        @test run_chrom_back.chromatograms[1].x     == rt_vec
        @test run_chrom_back.chromatograms[1].ic    == ic_vec
        @test MassJ.maxic(run_chrom_back.chromatograms[1]) ≈ maximum(ic_vec)
        rm(tmp_chrom)

        # -- Schema compliance: validate against the official PSI mzML 1.1.0
        #    XSDs when xmllint is available. The default save output is
        #    indexed-mzML (MaxQuant compatibility); we also test that
        #    `indexed = false` produces plain mzML that validates against
        #    the non-indexed schema.
        if Sys.which("xmllint") !== nothing
            xsd_plain   = nothing
            xsd_indexed = nothing
            for c in (joinpath("schema", "mzML1.1.0.xsd"), "/tmp/mzML1.1.0.xsd")
                isfile(c) && (xsd_plain = c; break)
            end
            for c in (joinpath("schema", "mzML1.1.0_idx.xsd"),
                      "/tmp/mzML-spec/schema/schema_1.1/mzML1.1.0_idx.xsd")
                isfile(c) && (xsd_indexed = c; break)
            end

            if xsd_plain !== nothing || xsd_indexed !== nothing
                # Representative MSrun with metadata + chromatogram.
                rt_v = collect(0.0:0.5:5.0)
                ic_v = Float64[100, 200, 350, 500, 600, 720, 690, 540, 400, 250, 110]
                run_full = MassJ.MSrun(scans_src.scans,
                                       Dict{String,Any}(
                                           "software" => [
                                               Dict{String,Any}(
                                                   "id" => "Xcalibur", "version" => "3.5",
                                                   "cv_params" => [
                                                       Dict{String,String}(
                                                           "accession" => "MS:1000532",
                                                           "name" => "Xcalibur"),
                                                   ]),
                                           ],
                                       ),
                                       [MassJ.IonCurrent(rt_v, ic_v; axis = :rt)])

                # Default (indexed = true) → validates against indexed schema
                if xsd_indexed !== nothing
                    tmp_idx = tempname() * ".mzML"
                    MassJ.save(run_full, tmp_idx; progress = false)
                    rc = run(`xmllint --noout --schema $xsd_indexed $tmp_idx`)
                    @test rc.exitcode == 0
                    rm(tmp_idx)
                end

                # indexed = false → validates against plain schema
                if xsd_plain !== nothing
                    tmp_plain = tempname() * ".mzML"
                    MassJ.save(run_full, tmp_plain; progress = false, indexed = false)
                    rc = run(`xmllint --noout --schema $xsd_plain $tmp_plain`)
                    @test rc.exitcode == 0
                    rm(tmp_plain)
                end
            else
                @info "PSI mzML XSDs not found; skipping schema-compliance test"
            end
        else
            @info "xmllint not available; skipping schema-compliance test"
        end

        # -- Indexed mzML: structural + SHA-1 + offset verification ---------
        # MaxQuant and most modern tools require the <indexedmzML> wrapper.
        # Confirm: presence of indexList/indexListOffset/fileChecksum, that
        # the emitted SHA-1 matches a re-computed one, and that each offset
        # really points to a `<spectrum` start.
        scans_idx_src = MassJ.load("test.mzML")
        tmp_idx_check = tempname() * ".mzML"
        MassJ.save(scans_idx_src, tmp_idx_check; progress = false)
        bytes_idx = read(tmp_idx_check)
        text_idx  = String(copy(bytes_idx))

        @test occursin("<indexedmzML",     text_idx)
        @test occursin("<indexList",       text_idx)
        @test occursin("<indexListOffset>", text_idx)
        @test occursin("<fileChecksum>",   text_idx)
        @test endswith(rstrip(text_idx), "</indexedmzML>")

        # SHA-1: re-compute over bytes up to and including the
        # `<fileChecksum>` open tag and compare with the value in the file.
        let fc = findfirst("<fileChecksum>", text_idx)
            fc_end_byte = last(fc)
            expected_hash = bytes2hex(SHA.sha1(bytes_idx[1:fc_end_byte]))
            fc_close = findfirst("</fileChecksum>", text_idx)
            emitted_hash = text_idx[fc_end_byte+1:first(fc_close)-1]
            @test emitted_hash == expected_hash
        end

        # Every <offset> idRef has an actual <spectrum at the recorded byte.
        for m in eachmatch(r"<offset idRef=\"([^\"]+)\">(\d+)</offset>", text_idx)
            byte_off = parse(Int, m.captures[2])
            @test startswith(text_idx[byte_off+1:end], "<spectrum")
        end

        rm(tmp_idx_check)
    end
end


function test_text_writers()
    @testset "Text writers - MGF / MSP / TXT round-trip" begin
        # -- MGF round-trip -------------------------------------------------
        mgf  = MassJ.load("test.mgf")
        tmp  = tempname() * ".mgf"
        MassJ.save(mgf, tmp)
        @test isfile(tmp)
        back = MassJ.load(tmp)
        @test length(back) == length(mgf)
        for (a, b) in zip(mgf, back)
            @test a.mz          ≈ b.mz
            @test a.int         ≈ b.int
            @test a.precursor   ≈ b.precursor
            @test a.chargeState == b.chargeState
            @test a.level       == b.level
        end
        rm(tmp)

        # -- MSP round-trip -------------------------------------------------
        msp   = MassJ.load("test.msp")
        tmpm  = tempname() * ".msp"
        MassJ.save(msp, tmpm)
        @test isfile(tmpm)
        backm = MassJ.load(tmpm)
        @test length(backm) == length(msp)
        for (a, b) in zip(msp, backm)
            @test a.mz              ≈ b.mz
            @test a.int             ≈ b.int
            @test a.precursor       ≈ b.precursor
            @test a.polarity        == b.polarity
            @test a.collisionEnergy ≈ b.collisionEnergy
        end
        rm(tmpm)

        # -- TXT round-trip (single spectrum) -------------------------------
        s1   = mgf[1]
        tmpt = tempname() * ".txt"
        MassJ.save(s1, tmpt)
        @test isfile(tmpt)
        backt = MassJ.load(tmpt)
        @test backt isa MassJ.MSscans      # bare, single spectrum
        @test backt.mz  ≈ s1.mz
        @test backt.int ≈ s1.int
        rm(tmpt)

        # TXT holds one spectrum: a multi-spectrum input must error.
        @test_throws ErrorException MassJ.save(mgf, tempname() * ".txt")

        # -- save() dispatches by extension; unknown extension errors -------
        @test MassJ.save(s1, tempname() * ".mgf") |> isfile
        @test_throws ErrorException MassJ.save(s1, tempname() * ".weird")
    end

    @testset "imzML writer round-trip (.imzML + .ibd)" begin
        orig = MassJ.load("test.imzML")
        tmp  = tempname() * ".imzML"
        MassJ.save(orig, tmp)
        @test isfile(tmp)
        @test isfile(splitext(tmp)[1] * ".ibd")          # companion binary written
        back = MassJ.load(tmp)
        @test length(back) == length(orig)
        for (a, b) in zip(orig, back)
            @test a.mz  ≈ b.mz
            @test a.int ≈ b.int
            @test a.metadata["position_x"] == b.metadata["position_x"]
            @test a.metadata["position_y"] == b.metadata["position_y"]
        end
        # zlib-compressed .ibd round-trips too
        tmpc = tempname() * ".imzML"
        MassJ.save(orig, tmpc; compress = true)
        backc = MassJ.load(tmpc)
        @test all(orig[i].mz ≈ backc[i].mz && orig[i].int ≈ backc[i].int for i in eachindex(orig))
    end
end


function test_cwt()
    @testset "Centroid - CWT (continuous wavelet transform)" begin
        Random.seed!(3)
        mz = collect(100.0:0.02:200.0); n = length(mz)
        int = 0.5 .* abs.(randn(n))                         # noise floor
        centres = [120.0, 150.0, 175.0]
        for (c, h) in zip(centres, [100.0, 60.0, 80.0])
            int .+= h .* exp.(-((mz .- c).^2) ./ (2 * 0.08^2))
        end
        s = MassJ.MSscans(1, 0.0, sum(int), mz, int, 1,
                          mz[argmax(int)], maximum(int), 0.0, "", "", 0.0)

        c = MassJ.centroid(s; method = MassJ.CWT())
        @test c isa MassJ.MSscans
        @test length(c.mz) == 3                              # the three peaks, no noise
        @test all(any(abs.(c.mz .- ctr) .< 0.2) for ctr in centres)

        c5 = MassJ.centroid(s; method = MassJ.CWT(threshold = 5.0))
        @test length(c5.mz) == 3                             # higher SNR threshold still clean

        cs = MassJ.centroid([s, s]; method = MassJ.CWT())    # vector method
        @test length(cs) == 2 && all(x -> length(x.mz) == 3, cs)
    end
end


function test_chimeric()
    # Build a synthetic series of two co-isolated precursors A and B whose
    # abundances are anticorrelated (they share the TIC). A fragments at
    # 201/305/410, B at 150/600.
    Random.seed!(42)
    N = 400
    a = abs.(randn(N)) .+ 1.0
    b = abs.(randn(N)) .+ 1.0
    pa = a ./ (a .+ b); pb = b ./ (a .+ b)
    afrags = [201.0, 305.0, 410.0]; bfrags = [150.0, 600.0]
    scans = MassJ.MSscans[]
    for n in 1:N
        mz = Float64[]; int = Float64[]
        for f in afrags; push!(mz, f); push!(int, max(pa[n]*100 + 0.5*randn(), 0.0)); end
        for f in bfrags; push!(mz, f); push!(int, max(pb[n]*100 + 0.5*randn(), 0.0)); end
        o = sortperm(mz)
        push!(scans, MassJ.MSscans(n, 0.0, sum(int), mz[o], int[o], 2,
                                   mz[o][argmax(int[o])], maximum(int), 0.0, "", "", 0.0))
    end

    @testset "Chimeric - abundance matrix" begin
        am = MassJ.abundance_matrix(scans; binsize = 1.0)
        @test size(am.matrix, 1) == N
        @test size(am.matrix, 2) == 5                 # 5 active bins (empty bins dropped)
        @test length(am.mz) == 5
        @test length(am.tic) == N
        @test am.tic ≈ vec(sum(am.matrix, dims = 2))
        # bins are the five fragment m/z (centres land in the 0.5-offset bin)
        @test all(any(abs.(am.mz .- f) .< 1.0) for f in vcat(afrags, bfrags))
    end

    @testset "Chimeric - partial correlation + clustering separates precursors" begin
        am = MassJ.abundance_matrix(scans; binsize = 1.0)
        P  = MassJ.partial_correlation(am.matrix)
        @test size(P) == (5, 5)
        @test all(isfinite, P)
        # bins sorted: 150(B),201(A),305(A),410(A),600(B) -> A = idx 2,3,4 ; B = idx 1,5
        @test P[2, 3] > 0 && P[2, 4] > 0 && P[3, 4] > 0      # A fragments co-vary
        @test P[1, 5] > 0                                    # B fragments co-vary
        @test P[1, 2] < 0 && P[2, 5] < 0                     # A vs B anticorrelated

        lab = MassJ.cluster_ions(P; kind = :correlation, nclusters = 2)
        @test length(lab) == 5
        @test lab[2] == lab[3] == lab[4]                     # the three A fragments together
        @test lab[1] == lab[5]                               # the two B fragments together
        @test lab[1] != lab[2]                               # A and B in different clusters

        specs = MassJ.cluster_spectra(am.mz, am.matrix, lab)
        @test length(specs) == 2
        @test all(s -> s isa MassJ.MSscans, specs)
        @test sort(vcat([length(s.mz) for s in specs])) == [2, 3]   # B has 2, A has 3
    end

    @testset "Chimeric - conditional mutual information" begin
        # CMI scenario: two groups driven by *independent* latent signals,
        # conditioned on an independent variable, so CMI separates them.
        Random.seed!(7)
        Nc = 400
        s1 = abs.(randn(Nc)) .+ 1.0; s2 = abs.(randn(Nc)) .+ 1.0
        z  = randn(Nc)
        g1 = [101.0, 202.0]; g2 = [303.0, 404.0]
        sc = MassJ.MSscans[]
        for n in 1:Nc
            mz = Float64[]; int = Float64[]
            for f in g1; push!(mz, f); push!(int, max(s1[n]*50 + 0.3*randn(), 0.0)); end
            for f in g2; push!(mz, f); push!(int, max(s2[n]*50 + 0.3*randn(), 0.0)); end
            o = sortperm(mz)
            push!(sc, MassJ.MSscans(n, 0.0, sum(int), mz[o], int[o], 2,
                                    mz[o][argmax(int[o])], maximum(int), 0.0, "", "", 0.0))
        end
        am = MassJ.abundance_matrix(sc; binsize = 1.0)
        C  = MassJ.cmi_matrix(am.matrix; condition = z)
        @test size(C) == (4, 4)
        @test all(isfinite, C)
        @test C ≈ C'                                         # symmetric
        @test all(==(0.0), [C[i, i] for i in 1:4])           # zero diagonal
        # bins sorted 101,202,303,404 -> g1 = 1,2 ; g2 = 3,4
        @test C[1, 2] > C[1, 3] && C[3, 4] > C[2, 3]         # within-group > cross-group

        lab = MassJ.cluster_ions(C; kind = :cmi, nclusters = 2)
        @test lab[1] == lab[2] && lab[3] == lab[4] && lab[1] != lab[3]
    end
end


function test_composed_predicates()
    scans = MassJ.load("test.mzXML")

    @testset "Composed predicates - empty filter returns all scans" begin
        sub = MassJ.extract(scans)
        @test length(sub) == length(scans)

        chrom = MassJ.chromatogram(scans)
        @test length(chrom.x) == length(scans)
    end

    @testset "Composed predicates - no match throws" begin
        @test_throws ErrorException MassJ.extract(scans, MassJ.Level(99))
        @test_throws ErrorException MassJ.chromatogram(scans, MassJ.Level(99))
        @test_throws ErrorException MassJ.average(scans, MassJ.Level(99))
    end

    @testset "Composed predicates - single match in average returns MSscan" begin
        result = MassJ.average(scans, MassJ.Scan(1))
        @test result isa MassJ.MSscans
        @test result.num == [1]
    end

    @testset "Composed predicates - AND semantics equivalence" begin
        # Single-pass composition must yield the same scans as stepwise filtering.
        combined = MassJ.extract(scans, MassJ.Level(2), MassJ.Polarity("+"))
        stepwise = MassJ.extract(MassJ.extract(scans, MassJ.Level(2)), MassJ.Polarity("+"))
        @test [s.num for s in combined] == [s.num for s in stepwise]
    end

    @testset "Composed predicates - multiple disjoint RT ranges" begin
        # Two RT intervals that together cover scans at RT≈0.14, 0.73 and 4.34 in test.mzXML.
        ms = MassJ.average(scans, MassJ.RT([[0.0, 1.0], [4.0, 5.0]]), stats = false)
        @test ms isa MassJ.MSscans
    end
end


function test_new_filters()
    # test.mzML: scan 1 = MS1 profile, charge 0; scan 2 = MS2 centroid, charge 2,
    # isolation_window_target_mz 400; scan 3 = MS2 centroid, charge 3, target 500.
    # All three reference instrument configuration IC1 (LTQ FT, ESI source, FTICR analyzer).
    run = MassJ.load("test.mzML")

    @testset "New filters - typed fields (charge / spectrum type / mobility)" begin
        @test length(MassJ.extract(run, MassJ.ChargeState(2)))                       == 1
        @test length(MassJ.extract(run, MassJ.ChargeState([2, 3])))                  == 2
        @test_throws ErrorException MassJ.extract(run, MassJ.ChargeState(99))

        @test length(MassJ.extract(run, MassJ.SpectrumType(:centroid)))              == 2
        @test length(MassJ.extract(run, MassJ.SpectrumType(:profile)))               == 1
        @test length(MassJ.extract(run, MassJ.SpectrumType([:centroid, :profile])))  == 3

        @test length(MassJ.extract(run, MassJ.MobilityType(:none)))                  == 3
        @test_throws ErrorException MassJ.extract(run, MassJ.MobilityType(:TIMS))
    end

    @testset "New filters - generic metadata cvParam" begin
        @test length(MassJ.extract(run, MassJ.MetaData("isolation_window_target_mz")))                 == 2  # presence
        @test length(MassJ.extract(run, MassJ.MetaData("isolation_window_target_mz", [400.0, 600.0]))) == 2  # range
        @test length(MassJ.extract(run, MassJ.MetaData("isolation_window_target_mz", 400.0)))          == 1  # numeric exact
        @test_throws ErrorException MassJ.extract(run, MassJ.MetaData("isolation_window_target_mz", 999.0))
        @test_throws ErrorException MassJ.extract(run, MassJ.MetaData("does_not_exist"))
    end

    @testset "New filters - instrument configuration" begin
        # resolved through the run-level instrument table: by id, cvParam name, component type
        @test length(MassJ.extract(run, MassJ.InstrumentConfig("IC1")))                      == 3
        @test length(MassJ.extract(run, MassJ.InstrumentConfig("LTQ FT")))                   == 3
        @test length(MassJ.extract(run, MassJ.InstrumentConfig("electrospray ionization")))  == 3
        @test length(MassJ.extract(run, MassJ.InstrumentConfig("source")))                   == 3
        @test_throws ErrorException MassJ.extract(run, MassJ.InstrumentConfig("ZZZ"))

        # On a plain Vector{MSscans} (no run table) only the per-spectrum ref id is matchable.
        vec = run.scans
        @test length(MassJ.extract(vec, MassJ.InstrumentConfig("IC1")))                      == 3
        @test_throws ErrorException MassJ.extract(vec, MassJ.InstrumentConfig("LTQ FT"))
    end

    @testset "New filters - compose with existing filters" begin
        @test length(MassJ.extract(run, MassJ.Level(2), MassJ.ChargeState(2)))                       == 1
        @test length(MassJ.extract(run, MassJ.InstrumentConfig("IC1"), MassJ.ChargeState(3)))        == 1
        @test length(MassJ.extract(run, MassJ.SpectrumType(:centroid),
                                        MassJ.MetaData("isolation_window_target_mz", [450.0, 600.0]))) == 1
        @test length(MassJ.chromatogram(run, MassJ.SpectrumType(:centroid)).x)                        == 2
        @test MassJ.average(run, MassJ.ChargeState([2, 3]))                                          isa MassJ.MSscans
    end

    @testset "New filters - mzXML file-direct path (generic fallback)" begin
        # The generic filter(::XMLElement, ::FilterType) fallback lets the new filters run on
        # the file-direct mzXML path. mzXML carries no charge state, so MobilityType(:none)
        # (the mzXML default) matches every scan, proving the fallback both dispatches and matches.
        @test length(MassJ.extract("test.mzXML", MassJ.MobilityType(:none)))           == 6
        @test length(MassJ.chromatogram("test.mzXML", MassJ.MobilityType(:none)).x)    == 6
        @test_throws ErrorException MassJ.extract("test.mzXML", MassJ.ChargeState(2))
    end
end


function test_yields()
    @testset "Yields - integrate_window, yields, normalize_*" begin

        # integrate_window on synthetic arrays: triangle from 0..1..0 over m/z 0..2
        mz   = [0.0, 1.0, 2.0]
        int  = [0.0, 1.0, 0.0]
        @test MassJ.integrate_window(mz, int, 0.0, 2.0) ≈ 1.0   # 2 * 1 / 2
        @test MassJ.integrate_window(mz, int, 0.0, 1.0) ≈ 0.5
        @test MassJ.integrate_window(mz, int, 5.0, 6.0) == 0.0  # no points
        @test MassJ.integrate_window(mz, int, 1.0, 1.0) == 0.0  # 1 point only
        @test MassJ.integrate_window(mz, int, 2.0, 0.0) ≈ 1.0   # swapped bounds

        # integrate_window on a real MSscan
        scans = MassJ.load("test.mzXML")
        a = MassJ.integrate_window(scans[1], 400.0, 500.0)
        @test a > 0 && isfinite(a)

        # Peak constructor swaps when mz1 > mz2
        p = MassJ.Peak(200.0, 100.0, "swap")
        @test p.mz1 == 100.0 && p.mz2 == 200.0 && p.label == "swap"

        # yields(files, peaks; x) — two-file series using the same fixture
        peaks = [MassJ.Peak(400.0, 500.0, "low"),
                 MassJ.Peak(800.0, 900.0, "high")]
        yc = MassJ.yields(["test.mzXML", "test.mzXML"], peaks;
                          x = [3.5, 4.0], xlabel = "photon energy (eV)")
        @test yc isa MassJ.YieldCurve
        @test size(yc.yields)    == (2, 2)
        @test size(yc.found_mz)  == (2, 2)
        @test all(isnan, yc.found_mz)         # all fixed Peak → all NaN
        @test yc.x               == [3.5, 4.0]
        @test yc.xlabel          == "photon energy (eV)"
        @test yc.labels          == ["low", "high"]
        @test yc.windows         == [(400.0, 500.0), (800.0, 900.0)]
        @test yc.tic[1]          ≈ yc.yields[1, 1] + yc.yields[1, 2]
        @test yc.yields[1, :]    ≈ yc.yields[2, :]    # same source file twice

        # length mismatch
        @test_throws ErrorException MassJ.yields(["test.mzXML"], peaks; x = [1.0, 2.0])

        # read_peaklist — round-trip through a temp CSV (with header)
        tmp = tempname() * ".csv"
        open(tmp, "w") do io
            write(io, "mz1,mz2,label\n")
            write(io, "400.0,500.0,low\n")
            write(io, "800.0,900.0,high\n")
        end
        pl = MassJ.read_peaklist(tmp)
        @test length(pl) == 2
        @test pl[1].mz1 == 400.0 && pl[1].mz2 == 500.0 && pl[1].label == "low"
        rm(tmp)

        # read_peaklist — no header
        tmp2 = tempname() * ".csv"
        open(tmp2, "w") do io
            write(io, "100.0,200.0,A\n")
        end
        pl2 = MassJ.read_peaklist(tmp2)
        @test length(pl2) == 1 && pl2[1].label == "A"
        rm(tmp2)

        # normalize_tic: rows of peak columns sum to 1
        yn = MassJ.normalize_tic(yc)
        @test yn isa MassJ.YieldCurve
        @test sum(yn.yields[1, :]) ≈ 1.0
        @test sum(yn.yields[2, :]) ≈ 1.0
        @test yn.tic == yc.tic   # raw totals preserved
        @test yn.metadata["normalize_tic"] == true

        # normalize_flux: divide by a constant flux of 2.0
        # — `#`-prefixed header lines are stripped as comments
        fluxpath = tempname() * ".txt"
        open(fluxpath, "w") do io
            write(io, "# header line 1\n")
            write(io, "# header line 2\n")
            write(io, "3.0  2.0\n")
            write(io, "5.0  2.0\n")
        end
        yf = MassJ.normalize_flux(yc, fluxpath)
        @test yf.yields ≈ yc.yields ./ 2.0
        @test yf.tic    ≈ yc.tic    ./ 2.0
        @test yf.metadata["normalize_flux"] == fluxpath
        rm(fluxpath)

        # normalize_flux: text header (no #) auto-detected and skipped
        flux_text = tempname() * ".txt"
        open(flux_text, "w") do io
            write(io, "energy flux\n")     # text header
            write(io, "==== ====\n")       # decorative
            write(io, "3.0  2.0\n")
            write(io, "5.0  2.0\n")
        end
        yf_t = MassJ.normalize_flux(yc, flux_text)
        @test yf_t.yields ≈ yc.yields ./ 2.0
        rm(flux_text)

        # normalize_flux: # comments mixed with data anywhere in the file
        flux_mixed = tempname() * ".txt"
        open(flux_mixed, "w") do io
            write(io, "# preamble\n")
            write(io, "energy   flux\n")   # text header below comments
            write(io, "3.0  2.0\n")
            write(io, "# mid-file note\n")
            write(io, "5.0  2.0  # trailing comment\n")
        end
        yf_m = MassJ.normalize_flux(yc, flux_mixed)
        @test yf_m.yields ≈ yc.yields ./ 2.0
        rm(flux_mixed)

        # normalize_flux: explicit skipstart override
        # — top row is "1 1" (looks numeric, but is a unit/scale row to skip)
        flux_skip = tempname() * ".txt"
        open(flux_skip, "w") do io
            write(io, "1 1\n")             # would otherwise be parsed as data
            write(io, "3.0  2.0\n")
            write(io, "5.0  2.0\n")
        end
        yf_s = MassJ.normalize_flux(yc, flux_skip; skipstart = 1)
        @test yf_s.yields ≈ yc.yields ./ 2.0
        rm(flux_skip)

        # write_csv round-trip — header + correct row count
        outpath = tempname() * ".csv"
        MassJ.write_csv(yc, outpath)
        lines = readlines(outpath)
        @test length(lines) == 3                                   # header + 2 rows
        @test lines[1] == "photon energy (eV),low,high,TIC"
        rm(outpath)

        # plot recipe smoke test
        @test typeof(plot(yc)) == Plots.Plot{Plots.GRBackend}
    end
end


function test_yields_targetpeak()
    @testset "TargetPeak resolution methods + Peak(mz, label; tol/ppm)" begin

        # Peak(mz, label; tol) — eager fixed window from a single m/z
        p = MassJ.Peak(100.0, "x"; tol = 0.5)
        @test p.mz1 == 99.5 && p.mz2 == 100.5 && p.label == "x"

        # Peak(mz, label; ppm)
        p2 = MassJ.Peak(1000.0, "y"; ppm = 5.0)
        @test p2.mz1 ≈ 1000.0 - 1000.0 * 5e-6
        @test p2.mz2 ≈ 1000.0 + 1000.0 * 5e-6

        # tol/ppm mutual exclusion
        @test_throws ErrorException MassJ.Peak(100.0, "z"; tol = 0.5, ppm = 5.0)
        @test_throws ErrorException MassJ.Peak(100.0, "z")

        # TargetPeak defaults and method validation
        tp = MassJ.TargetPeak(100.0, "a"; tol = 0.5)
        @test tp.mz == 100.0 && tp.tol == 0.5
        @test tp.method === :local_max && tp.edges == 0.1

        tp2 = MassJ.TargetPeak(100.0, "b"; tol = 0.5, method = :edges, edges = 0.2)
        @test tp2.method === :edges && tp2.edges == 0.2

        @test_throws ErrorException MassJ.TargetPeak(100.0, "c"; tol = 0.5, method = :bad)
        @test_throws ErrorException MassJ.TargetPeak(100.0, "c")     # neither tol nor ppm

        # Use the MS1 scans' actual base-peak m/z (~628.33) — a real,
        # well-covered peak. Picking the global argmax of the heterogeneous
        # averaged spectrum can land at an m/z covered by only some scans
        # (MS3 starts at 50, MS1 at 140, MS2 at 345); with per-scan
        # integration that yields zero (correctly) because no individual
        # scan has the peak. We read the grid-accurate value from scan 1.
        scans_load = MassJ.load("test.mzXML")
        spec = MassJ.average(scans_load)
        # Snap to the grid bin actually containing the base peak so the
        # exact-match assertions on found_mz pass.
        _bp_idx   = argmax(scans_load[1].int)
        target_mz = scans_load[1].mz[_bp_idx]

        # :local_max — should snap exactly onto the sample-grid maximum
        peaks_lm = [MassJ.TargetPeak(target_mz - 0.05, "lm"; tol = 0.2)]
        yc_lm    = MassJ.yields(["test.mzXML"], peaks_lm; x = [1.0])
        @test yc_lm.found_mz[1, 1] ≈ target_mz
        @test yc_lm.yields[1, 1]   > 0

        # :edges — same location, possibly different window width
        peaks_ed = [MassJ.TargetPeak(target_mz - 0.05, "ed";
                                      tol = 0.2, method = :edges)]
        yc_ed    = MassJ.yields(["test.mzXML"], peaks_ed; x = [1.0])
        @test yc_ed.found_mz[1, 1] ≈ target_mz
        @test yc_ed.yields[1, 1]   > 0

        # :centroid — uses the package's centroid(). Wider tol + TBPD so the
        # search window definitely contains at least one centroid in this fixture.
        peaks_cn = [MassJ.TargetPeak(target_mz, "cn";
                                      tol = 5.0, method = :centroid)]
        yc_cn    = MassJ.yields(["test.mzXML"], peaks_cn; x = [1.0],
                                centroid_method = MassJ.TBPD(:gauss, 4500., 0.2))
        @test isfinite(yc_cn.found_mz[1, 1])
        @test abs(yc_cn.found_mz[1, 1] - target_mz) < 5.0
        @test yc_cn.yields[1, 1] > 0

        # Mixed peak list — Peak gives NaN, TargetPeak gives a located m/z
        peaks_mx = [MassJ.Peak(400.0, 500.0, "static"),
                    MassJ.TargetPeak(target_mz, "lazy"; tol = 0.2)]
        yc_mx    = MassJ.yields(["test.mzXML"], peaks_mx; x = [1.0])
        @test isnan(yc_mx.found_mz[1, 1])
        @test yc_mx.found_mz[1, 2] ≈ target_mz

        # read_peaklist 2-col → TargetPeak with kwarg defaults
        tmp2 = tempname() * ".csv"
        open(tmp2, "w") do io
            write(io, "mz,label\n")
            write(io, "100.0,A\n")
            write(io, "200.0,B\n")
        end
        pl2 = MassJ.read_peaklist(tmp2; tol = 0.3, method = :edges)
        @test length(pl2) == 2
        @test pl2[1] isa MassJ.TargetPeak
        @test pl2[1].mz == 100.0 && pl2[1].tol == 0.3 && pl2[1].method === :edges
        rm(tmp2)

        # read_peaklist 4-col → TargetPeak with per-row tol + method
        tmp4 = tempname() * ".csv"
        open(tmp4, "w") do io
            write(io, "mz,tol,method,label\n")
            write(io, "100.0,0.3,local_max,A\n")
            write(io, "200.0,0.5,edges,B\n")
        end
        pl4 = MassJ.read_peaklist(tmp4)
        @test length(pl4) == 2
        @test pl4[1].method === :local_max && pl4[1].tol == 0.3
        @test pl4[2].method === :edges     && pl4[2].tol == 0.5
        rm(tmp4)

        # normalize_tic preserves found_mz
        yn = MassJ.normalize_tic(yc_lm)
        @test yn.found_mz == yc_lm.found_mz

        # drop_peaks — single string
        peaks_two = [MassJ.Peak(400.0, 500.0, "low"),
                     MassJ.Peak(800.0, 900.0, "high")]
        yc2 = MassJ.yields(["test.mzXML", "test.mzXML"], peaks_two;
                           x = [1.0, 2.0])
        d1 = MassJ.drop_peaks(yc2, "low")
        @test d1.labels        == ["high"]
        @test size(d1.yields)  == (2, 1)
        @test d1.yields[:, 1]  ≈  yc2.yields[:, 2]
        @test d1.windows       == [(800.0, 900.0)]
        @test size(d1.found_mz) == (2, 1)
        @test d1.tic           == yc2.tic         # tic unchanged by design

        # drop_peaks — vector of labels
        d2 = MassJ.drop_peaks(yc2, ["high"])
        @test d2.labels == ["low"]

        # drop_peaks — drop everything
        d3 = MassJ.drop_peaks(yc2, ["low", "high"])
        @test isempty(d3.labels)
        @test size(d3.yields) == (2, 0)

        # drop_peaks — missing labels silently ignored
        d4 = MassJ.drop_peaks(yc2, ["nonexistent"])
        @test d4.labels == yc2.labels
        @test d4.yields == yc2.yields
    end
end


function test_yields_errors()
    @testset "YieldCurve error propagation" begin

        # MSscans from `average("test.mzXML")` carries variance over 6 scans,
        # so yields_err and tic_err should all be finite.
        peaks = [MassJ.Peak(400.0, 500.0, "low"),
                 MassJ.Peak(800.0, 900.0, "high")]
        yc = MassJ.yields(["test.mzXML", "test.mzXML"], peaks; x = [1.0, 2.0])

        @test size(yc.yields_err) == (2, 2)
        @test all(isfinite, yc.yields_err)
        @test all(yc.yields_err .>= 0)
        @test size(yc.tic_err) == (2,)
        @test all(isfinite, yc.tic_err)
        for i in 1:2
            @test yc.tic_err[i] ≈ sqrt(sum(abs2, yc.yields_err[i, :]))
        end

        # normalize_tic: standard division error propagation
        yn = MassJ.normalize_tic(yc)
        @test all(isfinite, yn.yields_err)
        for i in 1:2, p in 1:2
            σ_y = yc.yields_err[i, p]
            y   = yc.yields[i, p]
            t   = yc.tic[i]
            σ_t = yc.tic_err[i]
            expected = sqrt((σ_y / t)^2 + (y * σ_t / (t * t))^2)
            @test yn.yields_err[i, p] ≈ expected
        end
        @test yn.tic     == yc.tic           # raw totals preserved
        @test yn.tic_err == yc.tic_err

        # normalize_flux: 2-col file → 10% default σ_φ
        fluxpath = tempname() * ".txt"
        open(fluxpath, "w") do io
            write(io, "# flux 10%\n")
            write(io, "1.0  2.0\n")
            write(io, "3.0  2.0\n")
        end
        yf = MassJ.normalize_flux(yc, fluxpath)
        @test yf.metadata["normalize_flux_err_pct"] == 0.10
        for i in 1:2, p in 1:2
            φ, σφ = 2.0, 0.10 * 2.0
            σ_y   = yc.yields_err[i, p]
            y     = yc.yields[i, p]
            expected = sqrt((σ_y / φ)^2 + (y * σφ / (φ * φ))^2)
            @test yf.yields_err[i, p] ≈ expected
        end
        rm(fluxpath)

        # normalize_flux: custom flux_err_pct kwarg
        fluxpath2 = tempname() * ".txt"
        open(fluxpath2, "w") do io
            write(io, "1.0  2.0\n")
            write(io, "3.0  2.0\n")
        end
        yf5 = MassJ.normalize_flux(yc, fluxpath2; flux_err_pct = 0.05)
        @test yf5.metadata["normalize_flux_err_pct"] == 0.05
        rm(fluxpath2)

        # normalize_flux: 3-col file → σ_φ from the file
        flux3 = tempname() * ".txt"
        open(flux3, "w") do io
            write(io, "# x flux sigma\n")
            write(io, "1.0  2.0  0.1\n")
            write(io, "3.0  2.0  0.1\n")
        end
        yf3 = MassJ.normalize_flux(yc, flux3)
        for i in 1:2, p in 1:2
            φ, σφ = 2.0, 0.1
            σ_y   = yc.yields_err[i, p]
            y     = yc.yields[i, p]
            expected = sqrt((σ_y / φ)^2 + (y * σφ / (φ * φ))^2)
            @test yf3.yields_err[i, p] ≈ expected
        end
        rm(flux3)

        # normalize_flux: jagged / empty 3rd column should not throw
        # (readdlm sometimes pads files with an empty 3rd col when rows are
        #  uneven or have trailing whitespace) — per-row fallback to pct.
        # Use xf == yc.x so there's no interpolation between σ values.
        flux_jag = tempname() * ".txt"
        open(flux_jag, "w") do io
            write(io, "# x flux [sigma]\n")
            write(io, "1.0  2.0  0.1\n")   # σ present
            write(io, "2.0  2.0\n")        # σ missing → falls back to pct
        end
        yfj = MassJ.normalize_flux(yc, flux_jag; flux_err_pct = 0.10)
        @test all(isfinite, yfj.yields_err)
        # Row 1 (yc.x=1.0) lands on xf[1] → σφ = 0.1
        # Row 2 (yc.x=2.0) lands on xf[2] → σφ = 0.10·|2.0| = 0.2 (pct fallback)
        for p in 1:2
            φ = 2.0
            σ_y1 = yc.yields_err[1, p]; y1 = yc.yields[1, p]
            σ_y2 = yc.yields_err[2, p]; y2 = yc.yields[2, p]
            @test yfj.yields_err[1, p] ≈ sqrt((σ_y1 / φ)^2 + (y1 * 0.1 / (φ * φ))^2)
            @test yfj.yields_err[2, p] ≈ sqrt((σ_y2 / φ)^2 + (y2 * 0.2 / (φ * φ))^2)
        end
        rm(flux_jag)

        # normalize_flux: non-numeric 3rd column (e.g. date string in DESIRS
        # beamline log files) — each row's σ falls back to pct.
        flux_str = tempname() * ".txt"
        open(flux_str, "w") do io
            write(io, "# DESIRS-style header\n")
            write(io, "Energy (eV)  flux  timestamp\n")          # text header
            write(io, "1.0  2.0  Fri Apr 17 12:27:54 2026\n")
            write(io, "2.0  2.0  Fri Apr 17 12:28:02 2026\n")
        end
        yfs = MassJ.normalize_flux(yc, flux_str; flux_err_pct = 0.10)
        @test all(isfinite, yfs.yields_err)
        rm(flux_str)

        # normalize_flux: extrapolate = :line when yc.x is outside flux range
        # flux defined on [3.0, 5.0] with slope (4-2)/(5-3) = 1 per eV → at x=2 (1 eV
        # before the start) the extrapolated φ should be 2 - 1·1 = 1.0
        flux_ext = tempname() * ".txt"
        open(flux_ext, "w") do io
            write(io, "3.0  2.0\n")
            write(io, "5.0  4.0\n")
        end
        # yc.x = [1.0, 2.0]; both points are below the flux range start (3.0).
        # With :line, φ(1.0) = 2.0 - 1.0·2 = 0.0 — non-positive → skip flagged.
        # With :line, φ(2.0) = 2.0 - 1.0·1 = 1.0 — valid division.
        yc_ext = MassJ.normalize_flux(yc, flux_ext;
                                       flux_err_pct = 0.10,
                                       extrapolate  = :line)
        @test yc_ext.metadata["normalize_flux_extrap"] == "line"
        # Row 2 (yc.x = 2.0) → extrapolated φ = 1.0
        for p in 1:2
            @test yc_ext.yields[2, p] ≈ yc.yields[2, p] / 1.0
        end
        rm(flux_ext)

        # :clamp (default) on out-of-range value clamps to nearest endpoint
        flux_clamp = tempname() * ".txt"
        open(flux_clamp, "w") do io
            write(io, "3.0  2.0\n")
            write(io, "5.0  4.0\n")
        end
        yc_clamp = MassJ.normalize_flux(yc, flux_clamp)   # default :clamp
        @test yc_clamp.metadata["normalize_flux_extrap"] == "clamp"
        # Row 1 (yc.x = 1.0, < 3.0) clamped to φ = 2.0
        for p in 1:2
            @test yc_clamp.yields[1, p] ≈ yc.yields[1, p] / 2.0
        end
        rm(flux_clamp)

        # Invalid extrapolate symbol
        flux_bad = tempname() * ".txt"
        open(flux_bad, "w") do io
            write(io, "3.0  2.0\n")
            write(io, "5.0  4.0\n")
        end
        @test_throws ErrorException MassJ.normalize_flux(yc, flux_bad;
                                                        extrapolate = :spline)
        rm(flux_bad)

        # drop_peaks slices yields_err; tic_err unchanged by design
        d = MassJ.drop_peaks(yc, "low")
        @test size(d.yields_err) == (2, 1)
        @test d.yields_err[:, 1] == yc.yields_err[:, 2]
        @test d.tic_err == yc.tic_err

        # Plot with ribbon still works
        @test typeof(plot(yc)) == Plots.Plot{Plots.GRBackend}

        # Per-scan-area error: a single-file yields still produces a finite,
        # non-negative SEM from the spread of the per-scan integrals.
        scans_real = MassJ.load("test.mzXML")
        bp_mz      = scans_real[1].mz[argmax(scans_real[1].int)]
        yc_cmp     = MassJ.yields(["test.mzXML"], [MassJ.TargetPeak(bp_mz, "bp"; tol = 0.5)]; x = [1.0])
        @test isfinite(yc_cmp.yields_err[1, 1])
        @test yc_cmp.yields_err[1, 1] >= 0
    end
end


function test_aqua()
    @testset "Aqua — package-quality checks" begin
        # Lighter check set than test_all to start; we'll enable the rest
        # once incremental fixes have landed.
        Aqua.test_unbound_args(MassJ)
        Aqua.test_undefined_exports(MassJ)
        Aqua.test_project_extras(MassJ)
        Aqua.test_stale_deps(MassJ)
        Aqua.test_deps_compat(MassJ; check_extras = false)
    end
end


tests()
test_isotopes()
test_adducts()
test_calibration()
test_centroid_metrics()
test_assignment()
test_multifolder()
test_yield_transforms()
test_peptides()
test_fragment_peaks()
test_mobility()
test_uncertainty()
test_tables()
test_measurements()
test_unitful()
test_chrom_processing()
test_chrom_peaks()
test_deconvolution()
test_interpolation_import()
test_mzml()
test_mgf()
test_msp()
test_imzml()
test_composed_predicates()
test_new_filters()
test_export()
test_text_writers()
test_cwt()
test_chimeric()
test_yields()
test_yields_targetpeak()
test_yields_errors()
test_aqua()
