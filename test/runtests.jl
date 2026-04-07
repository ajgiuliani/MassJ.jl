using MSj, Test
using Plots
using DataStructures

function tests()
    @testset "Subset of tests"  begin
        inf = MSj.info("test.mzXML", verbose = true)
        @test inf[1] == "parentFile: test.raw"                                         #1
        @test inf[9] == "6 scans"                                                      #2
        @test inf[10] == "MS1+"                                                        #3
        @test inf[11] == "MS2+ 1255.5  CID(CE=18)"                                     #4
        @test inf[12] == "MS3+ 902.33  PQD(CE=35)"                                     #5

        scans = MSj.load("test.mzXML")
        @test eltype(scans)              == MSj.MSscan                                 #6
        @test length(scans)              == 6                                          #7
        @test scans[1].num               == 1                                          #8
        @test scans[2].level             == 2                                          #9
        @test scans[3].polarity          == "+"                                        #10
        @test scans[2].activationMethod  == "CID"                                      #11
        @test scans[3].collisionEnergy   == 35.0                                       #12
        @test size(scans[1].int, 1)      == 22320                                      #13

        rt = MSj.retention_time("test.mzXML")
        @test length(rt) == 6                                                          #14

        cr = MSj.chromatogram("test.mzXML", method = MSj.TIC() )
        @test length(cr.rt) == 6                                                       #15

        cr = MSj.chromatogram("test.mzXML", method = MSj.MZ([0, 500]))
        @test length(cr.rt) == 6                                                       #16

        cr = MSj.chromatogram("test.mzXML", method = MSj.∆MZ([1000, 1]))
        @test length(cr.rt) == 6                                                       #17

        cr = MSj.chromatogram("test.mzXML", method = MSj.BasePeak())
        @test length(cr.rt) == 6                                                       #18

        cr = MSj.chromatogram("test.mzXML", MSj.Polarity("+"), MSj.Scan(2),MSj.Precursor(1255.5), MSj.Activation_Energy(18), MSj.Activation_Method("CID"), MSj.Level(2) )
        @test length(cr.rt) == 1                                                       #19

        rt = MSj.retention_time(scans)
        @test length(rt) == 6                                                          #20

        cr = MSj.chromatogram(scans, method = MSj.TIC() )
        @test length(cr.rt) == 6                                                       #21

        cr = MSj.chromatogram(scans, method = MSj.MZ([0, 500]))
        @test length(cr.rt) == 6                                                       #22

        cr = MSj.chromatogram(scans, method = MSj.∆MZ([1000, 1]))
        @test length(cr.rt) == 6                                                       #23

        cr = MSj.chromatogram(scans, method = MSj.BasePeak())
        @test length(cr.rt) == 6                                                       #24

        cr = MSj.chromatogram(scans, MSj.Polarity("+"),MSj.Scan(2),MSj.Precursor(1255.5),MSj.Activation_Energy(18),MSj.Activation_Method("CID"),MSj.Level(2) )
        @test (cr.rt, cr.ic) == ([0.7307], [9727.2])                                   #25

        cr = MSj.chromatogram(scans, MSj.Polarity(["+"]),MSj.Scan([2,3]),MSj.Precursor([1255.5, 902.33]),MSj.Activation_Energy([18, 35]),MSj.Activation_Method(["CID", "PQD"]),MSj.Level([2, 3]) )
        @test (cr.rt, cr.ic) == ([0.7307, 2.1379], [9727.2, 11.3032])                  #26

        ms = MSj.average("test.mzXML")
        @test length(ms.num) == 6                                                      #27

        ms = MSj.average("test.mzXML", MSj.Polarity("+"),MSj.Scan(2),MSj.Precursor(1255.5),MSj.Activation_Energy(18),MSj.Activation_Method("CID"),MSj.RT(1),MSj.IC([0, 1e4]))
        @test ms isa MSj.MSscan                                                        #28
        @test ms.num == 2                                                              #29

        ms = MSj.average(scans)
        @test length(ms.num) == 6                                                      #30

        ms = MSj.average(scans, MSj.Polarity("+"),MSj.Scan(2),MSj.Precursor(1255.5),MSj.Activation_Energy(18),MSj.Activation_Method("CID"),MSj.RT(1),MSj.IC([0, 1e4]))
        @test ms isa MSj.MSscan                                                        #31
        @test ms.num == 2                                                              #32

        ms = MSj.average(scans, MSj.Polarity(["+"]),MSj.Scan([2,3]),MSj.Precursor([1255.5, 902.33]),MSj.Activation_Energy([18, 35]),MSj.Activation_Method(["CID", "PQD"]),MSj.RT([1,2]),MSj.IC([0, 1e4]))
        @test ms isa MSj.MSscans                                                       #33
        @test ms.num == [2, 3]                                                         #34

        ms = MSj.average("test.mzXML", MSj.RT( [[1,2], [2,3]] ), stats = false )
        @test ms isa MSj.MSscans                                                       #35
        @test ms.num == [2, 3, 4]                                                      #36

        ms = MSj.average("test.mzXML", MSj.Polarity(["+"]),MSj.Scan([2,3]),MSj.Precursor([1255.5, 902.33]),MSj.Activation_Method(["CID", "PQD"]),MSj.RT([1,2]),MSj.IC([0, 1e4]))   #MSj.Activation_Energy([18., 35.]),
        @test ms isa MSj.MSscans                                                       #37
        @test ms.num == [2, 3]                                                         #38

        cr = MSj.chromatogram("test.mzXML", MSj.Polarity(["+"]),MSj.Scan([2,3]),MSj.Precursor([1255.5, 902.33]),MSj.Activation_Method(["CID", "PQD"]),MSj.Level([2, 3]) )   #MSj.Activation_Energy([18.0, 35.0]),
        @test length(cr.rt) == 2                                                       #39

        ms = MSj.average(scans, MSj.RT( [[1,2], [2,3]] ), stats = false )
        @test ms isa MSj.MSscans                                                       #40
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

        a = MSj.avg(scans[1], scans[2])
        @test a.num == [1,2]                                                           #55

        a = MSj.avg(scans[1], scans[4])
        @test a.num == [1,4]                                                           #56

        info = MSj.info("test64.mzXML")
        @test info[2] == "MS1-"                                                        #57

        scans = MSj.load("test64.mzXML")
        @test eltype(scans)              == MSj.MSscan                                 #58

        info = MSj.info("test.mzXMLM")
        @test info.msg == "File format not supported."                                 #59

        scans = MSj.load("test.mzXMLL")
        @test info.msg == "File format not supported."                                 #60

        scans = MSj.load("bad1.mzXML")
        @test scans.msg == "Not an mzXML file."                                        #61

        scans = MSj.info("bad1.mzXML")
        @test scans.msg == "Not an mzXML file."                                        #62

        scans = MSj.load("bad2.mzXML")
        @test scans[1].num == 0                                                        #63

        scans = MSj.load("bad3.mzXML")
        @test scans[1].num == scans[2].num == 0                                        #64

        cr = MSj.chromatogram("test.mzXML", method = MSj.∆MZ([1, 2]))
        @test cr.msg == "Bad mz ± ∆mz values."                                         #65

        cr = MSj.chromatogram(scans, method = MSj.∆MZ([1, 2]))
        @test cr.msg == "Bad mz ± ∆mz values."                                         #66

        scans = MSj.load("test.mzXML")
        @test MSj.smooth(scans[1], method = MSj.SG(7,15,0)) isa MSj.MSscan             #67

        a = MSj.avg(scans[1], scans[4])
        @test MSj.smooth(a) isa MSj.MSscans                                            #68

        a = scans[1] * scans[4]
        @test a.num == [1,4]                                                           #69

        a = (scans[2]+scans[3]) - scans[1]
        @test a.num == [2, 3]                                                          #70

        a = (scans[1] + scans[4]) - (scans[1] - scans[4])
        @test a.num == [1, 4]                                                          #71

        a = scans[1] + MSj.avg(scans[2], scans[5])
        @test a.num == [1, 2, 5]                                                       #72

        a = MSj.smooth(scans[1], method = MSj.SG(5,9,0))
        @test a.num == 1                                                               #73

       a = MSj.centroid(scans[1], method = MSj.TBPD(:gauss, 4500., 0.2))               #74
       @test length(a.int) == 957

       @test typeof(plot(scans[1], method = :relative)) == Plots.Plot{Plots.GRBackend} #75
       @test typeof(plot(scans[1], method = :absolute)) == Plots.Plot{Plots.GRBackend} #76

       a = MSj.avg(scans[2], scans[5])
       @test typeof(plot( a, method = :relative )) == Plots.Plot{Plots.GRBackend}      #77
       @test typeof(plot( a, method = :absolute )) == Plots.Plot{Plots.GRBackend}      #78

       cr = MSj.chromatogram(scans)
       @test typeof(plot( cr, method = :relative )) == Plots.Plot{Plots.GRBackend}     #79
       @test typeof(plot( cr, method = :absolute )) == Plots.Plot{Plots.GRBackend}     #80

       a = MSj.centroid(scans[1], method = MSj.TBPD(:voigt, 4500., 0.2))
       @test length(a.int) == 961                                                      #81

       a = MSj.centroid(scans[1], method = MSj.TBPD(:lorentz, 4500., 0.2))
       @test length(a.int) == 964                                                      #82

       a = MSj.centroid(scans[1], method = MSj.TBPD(:other, 4500., 0.2))
       @test a.msg == "Unsupported peak profile. Use :gauss, :lorentz or :voigt."      #83

       a = MSj.centroid(scans[1], method = MSj.SNRA(1., 100))
       @test length(a.int) == 109                                                      #84

       s1 = MSj.extract(scans, MSj.Activation_Energy([18,35]))
       @test length(s1) == 4                                                           #85

       s1 = MSj.extract("test.mzXML", MSj.Activation_Energy(18))
       @test length(s1) == 2                                                           #86

       s1 = MSj.extract(scans, MSj.Scan(1))
       @test length(s1) == 1                                                           #87

       s1 = MSj.extract("test.mzXML", MSj.Scan(1))
       @test length(s1) == 1                                                           #88

       bs = MSj.baseline_correction(scans, method = MSj.TopHat(1))
       @test length(bs) == 6                                                           #89

       bs = MSj.baseline_correction(scans[1], method = MSj.TopHat(1))
       @test length(bs.int) == length(scans[1].int)                                    #90

       c = MSj.centroid(scans, method = MSj.TBPD(:gauss, 4500., 0.2)) ;
       bs = MSj.baseline_correction( c, method = MSj.LOESS(3))
       @test length(bs) == 6                                                           #91

       bs = MSj.baseline_correction(c[1], method = MSj.LOESS(3))
       @test length(bs.int) == length(c[1].int)                                        #92

       bs = MSj.baseline_correction(scans, method = MSj.IPSA(51,100))
       @test length(bs) == 6                                                           #93

       bs = MSj.baseline_correction(scans[1], method = MSj.IPSA(51,100))
       @test length(bs.int) == length(scans[1].int)                                    #94

       a = MSj.smooth(scans, method = MSj.SG(5,9,0))
       @test length(a) == 6                                                            #95

       c = MSj.centroid(scans, method = MSj.TBPD(:lorentz, 4500., 0.2)) ;
       d = MSj.centroid(scans, method = MSj.TBPD(:voigt, 4500., 0.2)) ;
       @test length(c) == length(d)                                                    #96

       a = MSj.centroid(scans[3], method = MSj.SNRA(1., 100))
       @test length(a.int) == 0                                                        #97

       bs = MSj.baseline_correction(scans[1], method = MSj.IPSA(50,100))
       @test length(bs.int) == length(scans[1].int)                                    #98

       cr = MSj.chromatogram(scans, method = MSj.BasePeak() )
       @test length(cr.rt) == 6                                                        #99

       f = MSj.formula("CH3(13C)10H3Kr(NaH2)2")                                        #100
      @test f == Dict("Na" => 2,"Kr" => 1,"C" => 1,"13C" => 10,"H" => 10)

       m = MSj.masses("C254 H377 N65 O75 S6")                                          #101
      @test m == Dict("Monoisotopic" => 5729.60087099839, "Average" => 5733.55, "Nominal" => 5727.0)

      I = MSj.isotopic_distribution("CH4", 0.9999, charge = +1)                        #102
      @test I[2,1:end] == [16.03130012908, 0.9887541751052761, 1, 0, 4, 0]

      a = MSj.simulate(I, 0.4, Npoints = 5)                                            #103
      @test a.int == [100.0, 6.035851011021856, 0.06994625998243831, 1.1368602493290418, 0.06861901924947518]

       m = MSj.masses(f)                                                               #104
      @test m == Dict("Monoisotopic" => 282.0028349717, "Average" => 281.902086912, "Nominal" => 282.0)

      I = MSj.isotopic_distribution(f, 0.9999, charge = +1)                            #105
      @test I[2,1:end][1:2] == [282.0028349717, 0.5630635281692917]

      a = MSj.simulate(I, 0.4, model=:lorentz, Npoints = 5)                            #106
      @test a.int == [0.5359868694750152, 16.332387108915455, 100.0, 0.560482304632663, 0.0834376623225204]

    end
end


function test_isotopes()
    @testset "Isotopes - type stability and optimizations" begin

        # Elements dict is const (enables type inference)
        @test isconst(MSj, :Elements)

        # PriorityQueues are typed
        pq = PriorityQueue{Vector{Int},Float64}()
        pq[[1,0]] = 0.5
        pq[[0,1]] = 0.3
        @test peek(pq) == ([0,1] => 0.3)   # min-heap: lowest first

        # formula parsing
        @test MSj.formula("C2H6O") == Dict("C" => 2, "H" => 6, "O" => 1)
        @test MSj.formula("H2O") == Dict("H" => 2, "O" => 1)
        @test MSj.formula("NaCl") == Dict("Na" => 1, "Cl" => 1)
        @test_throws ErrorException MSj.formula("123bad")
        @test_throws ErrorException MSj.formula("Xx")

        # masses calculation
        m = MSj.masses("H2O")
        @test m["Monoisotopic"] ≈ 18.01056468474
        @test m["Average"] ≈ 18.015
        @test m["Nominal"] == 18.0

        m2 = MSj.masses(Dict("H" => 2, "O" => 1))
        @test m2 == m

        # stirling approximation vs exact log factorial
        @test MSj.stirling(100) ≈ log(factorial(big(100))) atol=0.01
        @test MSj.stirling(500) > 0

        # isotopologue_probability - low mass path (Natoms < 20)
        prob_H2 = MSj.isotopologue_probability(Dict("H" => 2), Dict("H" => [2, 0]), MSj.Elements)
        @test prob_H2 ≈ MSj.Elements["H"][1].f^2 atol=1e-10

        # isotopologue_probability - high mass path (Natoms >= 20)
        prob_C20 = MSj.isotopologue_probability(Dict("C" => 20), Dict("C" => [20, 0]), MSj.Elements)
        @test prob_C20 ≈ 0.8049835604738165 atol=1e-8

        # isotopologue_mass
        cm = MSj.isotopologue_mass([Pair("H", 2)], Dict("H" => [2, 0]), MSj.Elements)
        @test cm ≈ 2 * MSj.Elements["H"][1].m

        # most_probable_isotopologue
        alpha = MSj.most_probable_isotopologue(Dict("C" => 10, "H" => 22), MSj.Elements)
        @test sum(alpha["C"]) == 10     # conservation
        @test sum(alpha["H"]) == 22     # conservation
        @test alpha["C"][1] >= alpha["C"][2]  # 12C more abundant than 13C

        # hill_climbing finds optimum
        f_obj(x) = -(x[1] - 3)^2 - (x[2] - 2)^2
        P = MSj.hill_climbing([1, 4], f_obj)
        @test P == [3, 2]

        # hill_climbing with single-element vector (early return)
        P_single = MSj.hill_climbing([5], x -> -x[1]^2)
        @test P_single == [5]  # no neighbors possible

        # isotopic distribution - basic
        I = MSj.isotopic_distribution("H2O", 0.99, charge = 1)
        @test size(I, 2) == 7   # Masses, Probability, + isotope columns
        @test I[2, 1] ≈ 18.01056468474   # monoisotopic mass
        @test I[2, 2] ≈ 0.9973367663173334  # probability
        @test I[1, 1] == "Masses"
        @test I[1, 2] == "Probability"

        # isotopic distribution - charge state divides mass
        I2 = MSj.isotopic_distribution("H2O", 0.99, charge = 2)
        @test I2[2, 1] ≈ 18.01056468474 / 2

        # isotopic distribution - larger molecule
        I3 = MSj.isotopic_distribution("C254 H377 N65 O75 S6", 0.5)
        @test size(I3, 1) > 2          # more than just header
        @test I3[1, 1] == "Masses"
        probs = [I3[i, 2] for i in 2:size(I3, 1)]
        @test all(p -> p > 0, probs)   # all probabilities positive
        @test sum(probs) ≈ 0.5 atol=0.1

    end
end


function test_deconvolution()
    @testset "Deconvolution - helpers and integration" begin

        # is_evenly_spaced
        @test MSj.is_evenly_spaced([1.0, 2.0, 3.0, 4.0]) == true
        @test MSj.is_evenly_spaced([1.0, 2.0, 3.5, 4.0]) == false
        @test MSj.is_evenly_spaced(collect(range(100.0, stop=200.0, length=100))) == true

        # resampling
        X = [1.0, 2.0, 3.5, 5.0]
        Y = [10.0, 20.0, 35.0, 50.0]
        newX, newY = MSj.resampling(X, Y)
        @test length(newX) == length(X)
        @test newX[1] ≈ X[1]
        @test newX[end] ≈ X[end]
        @test MSj.is_evenly_spaced(newX)

        # new_mass
        @test MSj.new_mass(500.0, 5, 1, 1.00782503227) ≈ 416.83463750537834
        @test MSj.new_mass(500.0, 5, -1, 1.00782503227) ≈ 624.7480437419325
        # new_mass with zero shift returns original
        @test MSj.new_mass(500.0, 5, 0, 1.0) ≈ 500.0

        # get_peak_shape
        mz = collect(range(100.0, stop=200.0, length=1000))
        g = MSj.get_peak_shape(MSj.gauss, 0.5, mz)
        @test length(g) > 0
        @test sum(g) ≈ 1.0 atol=1e-10      # normalized
        @test g[argmax(g)] > 0              # has a peak
        # approximately symmetric shape
        @test g[1] ≈ g[end] rtol=1.0

        # get_peak_shape with lorentz
        g_lor = MSj.get_peak_shape(MSj.lorentz, 0.5, mz)
        @test sum(g_lor) ≈ 1.0 atol=1e-10

        # figure_of_merit
        h = [1.0, 2.0, 3.0, 4.0, 5.0]
        @test MSj.figure_of_merit(h, h) ≈ 1.0           # perfect match
        @test MSj.figure_of_merit(h, h .+ 0.1) ≈ 0.995  # close match
        @test MSj.figure_of_merit(h, h .+ 0.1) > 0.99

        # chargefilter - basic shape and properties
        mz_test = collect(range(400.0, stop=600.0, length=200))
        int_test = zeros(200)
        int_test[100] = 100.0   # single peak at ~500
        f_test = zeros(10, 200)
        for i in 1:10
            f_test[i, :] = int_test
        end
        s = MSj.chargefilter(mz_test, int_test, f_test, (1, 10), 1.0, 1)
        @test size(s) == (10, 200)
        @test all(s .>= 0)     # non-negative output

        # project_N_convolve
        g_small = MSj.get_peak_shape(MSj.gauss, 2.0, mz_test)
        c = MSj.project_N_convolve(s, g_small)
        @test length(c) == 200
        @test all(c .>= 0)

        # Charges struct construction
        ch = MSj.Charges(adduct="H", range=(5, 15), width=2)
        @test ch.adduct == "H"
        @test ch.range == (5, 15)
        @test ch.width == 2

        # Charges with default width
        ch2 = MSj.Charges(adduct="Na", range=(1, 5))
        @test ch2.width == 1

        # _resolve_shape
        scans = MSj.load("test.mzXML")
        model_gauss = MSj._resolve_shape(scans[1], :gauss)
        @test model_gauss === MSj.gauss
        model_lor = MSj._resolve_shape(scans[1], :lorentz)
        @test model_lor === MSj.lorentz
        model_voigt = MSj._resolve_shape(scans[1], :voigt)
        @test model_voigt === MSj.voigt
        @test_throws ErrorException MSj._resolve_shape(scans[1], :invalid)

        # _resolve_FWHM with explicit values
        @test MSj._resolve_FWHM(scans[1], MSj.gauss, -1, 0.5) == 0.5   # explicit FWHM
        @test MSj._resolve_FWHM(scans[1], MSj.gauss, 1000.0, -1) == 0.5  # R=1000 → 500/1000=0.5
        @test_throws ErrorException MSj._resolve_FWHM(scans[1], MSj.gauss, 1000.0, 0.5)  # both

        # roughguess_FWHM runs without error
        w = MSj.roughguess_FWHM(scans[1])
        @test w > 0

        # deconv method exists with correct dispatch
        @test hasmethod(MSj.deconv, Tuple{MSj.MSscan, MSj.Charges})
        @test hasmethod(MSj.deconv, Tuple{MSj.MSscans, MSj.Charges})

    end
end


function test_interpolation_import()
    @testset "Interpolations import fix" begin
        # Line is accessible from MSj (needed for extrapolation)
        @test isdefined(MSj, :Line)
        @test isdefined(MSj, :LinearInterpolation)

        # Verify interpolation with Line extrapolation works
        x = [1.0, 2.0, 3.0, 4.0]
        y = [10.0, 20.0, 30.0, 40.0]
        itp = MSj.LinearInterpolation(x, y, extrapolation_bc=MSj.Line())
        @test itp(2.5) ≈ 25.0
        @test itp(0.0) ≈ 0.0    # extrapolated
    end
end


function test_mzml()
    @testset "mzML format" begin
        # info
        inf = MSj.info("test.mzML")
        @test any(contains(s, "3 scans") for s in inf)
        @test any(contains(s, "MS1+") for s in inf)
        @test any(contains(s, "MS2+") for s in inf)

        inf_v = MSj.info("test.mzML", verbose=true)
        @test any(contains(s, "test.raw") for s in inf_v)
        @test any(contains(s, "LTQ FT") for s in inf_v)

        # load all
        scans = MSj.load("test.mzML")
        @test length(scans) == 3
        @test eltype(scans) == MSj.MSscan

        # Scan 1: MS1+, profile, rt=0.5 min
        s1 = scans[1]
        @test s1.level == 1
        @test s1.rt ≈ 0.5
        @test s1.polarity == "+"
        @test s1.spectrumType == :profile
        @test s1.tic ≈ 19000.0
        @test s1.basePeakMz ≈ 400.0
        @test s1.basePeakIntensity ≈ 8000.0
        @test s1.precursor ≈ 0.0
        @test s1.chargeState == 0
        @test length(s1.mz) == 5
        @test s1.mz ≈ [100.0, 200.0, 300.0, 400.0, 500.0]
        @test s1.int ≈ [1000.0, 5000.0, 3000.0, 8000.0, 2000.0]

        # Scan 2: MS2+ CID, centroid, rt=1.0 min (60s converted)
        s2 = scans[2]
        @test s2.level == 2
        @test s2.rt ≈ 1.0
        @test s2.polarity == "+"
        @test s2.spectrumType == :centroid
        @test s2.precursor ≈ 400.0
        @test s2.chargeState == 2
        @test s2.activationMethod == "CID"
        @test s2.collisionEnergy ≈ 25.0
        @test length(s2.mz) == 4
        @test s2.mz ≈ [110.0, 150.0, 200.0, 250.0]

        # Scan 3: MS2+ HCD
        s3 = scans[3]
        @test s3.level == 2
        @test s3.rt ≈ 1.5
        @test s3.precursor ≈ 500.0
        @test s3.chargeState == 3
        @test s3.activationMethod == "HCD"
        @test s3.collisionEnergy ≈ 30.0
        @test length(s3.mz) == 3

        # retention_time
        rt = MSj.retention_time("test.mzML")
        @test length(rt) == 3
        @test rt ≈ [0.5, 1.0, 1.5]

        # chromatogram (TIC)
        chrom = MSj.chromatogram("test.mzML")
        @test length(chrom.rt) == 3
        @test chrom.rt ≈ [0.5, 1.0, 1.5]
        @test chrom.ic ≈ [19000.0, 4800.0, 2100.0]
        @test chrom.maxic ≈ 19000.0

        # chromatogram with filter
        chrom2 = MSj.chromatogram("test.mzML", MSj.Level(2))
        @test length(chrom2.rt) == 2

        # chromatogram base peak
        chrom_bp = MSj.chromatogram("test.mzML", method=MSj.BasePeak())
        @test chrom_bp.ic ≈ [8000.0, 2000.0, 1200.0]

        # extract
        ms2 = MSj.extract("test.mzML", MSj.Level(2))
        @test length(ms2) == 2
        @test ms2[1].precursor ≈ 400.0
        @test ms2[2].precursor ≈ 500.0

        # average
        avg = MSj.average("test.mzML", MSj.Level(2))
        @test avg isa MSj.MSscans

        # New fields present
        @test s1.mobilityType == :none
        @test s1.driftTime ≈ -1.0
        @test s1.compensationVoltage ≈ 0.0
        @test s1.metadata isa Dict{String,Any}
    end
end


function test_mgf()
    @testset "MGF format" begin
        # info
        inf = MSj.info("test.mgf")
        @test inf[1] == "3 scans"
        @test any(contains(s, "400.0") for s in inf)
        @test any(contains(s, "500.0") for s in inf)
        @test any(contains(s, "600.0") for s in inf)

        # load all
        scans = MSj.load("test.mgf")
        @test length(scans) == 3
        @test eltype(scans) == MSj.MSscan

        # Scan 1
        s1 = scans[1]
        @test s1.num == 1        # sequential index
        @test s1.rt ≈ 0.5       # 30s / 60
        @test s1.level == 2      # MGF default
        @test s1.precursor ≈ 400.0
        @test s1.chargeState == 2
        @test s1.polarity == "+"
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
        @test s2.num == 2
        @test s2.rt ≈ 1.0
        @test s2.precursor ≈ 500.0
        @test s2.chargeState == 3

        # Scan 3
        s3 = scans[3]
        @test s3.num == 3
        @test s3.precursor ≈ 600.0
        @test s3.chargeState == 1
        @test length(s3.mz) == 4

        # retention_time
        rt = MSj.retention_time("test.mgf")
        @test length(rt) == 3
        @test rt ≈ [0.5, 1.0, 1.5]

        # chromatogram
        chrom = MSj.chromatogram("test.mgf")
        @test length(chrom.rt) == 3
        @test chrom.maxic ≈ 4800.0

        # extract (filter by precursor)
        sub = MSj.extract("test.mgf", MSj.Precursor(500.0))
        @test length(sub) == 1
        @test sub[1].precursor ≈ 500.0

        # average
        avg = MSj.average("test.mgf")
        @test avg isa MSj.MSscans

        # New MSscan fields
        @test s1.mobilityType == :none
        @test s1.driftTime ≈ -1.0
        @test s1.compensationVoltage ≈ 0.0
    end
end


tests()
test_isotopes()
test_deconvolution()
test_interpolation_import()
test_mzml()
test_mgf()
