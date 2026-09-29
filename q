[33mcommit 7677488e086659e779f11e8c90ceadca6cf48626[m[33m ([m[1;36mHEAD -> [m[1;32mmaster[m[33m, [m[1;31morigin/master[m[33m, [m[1;31morigin/HEAD[m[33m)[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Tue Jun 15 15:55:57 2021 -0700

    Verify Fresnel propagation
    
    Add <Debug_Fresnel_Propagation.ipynb>.  Add `verbose` option in
    sleec_input.import_pars().  Modified <sleec_alter._test_dm_disk_harmonic()>.

[33mcommit feb59fbf91a10189de32dc1454b7b35b8eae2195[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Wed Jun 9 14:18:00 2021 -0700

    Add function to reset DM surface(s) using disk harmonic basis set
    
    New code passes <_test_dm_disk_harmonic()>.  However, notice that phase
    plots at DM2 and relay DM always show flat values; need to investigate.

[33mcommit 31ed803e89e52ff798c8eebac0d19e57058db038[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Tue Jun 8 15:42:55 2021 -0700

    Correct function doc string

[33mcommit 807148b7651b6fa7556070872a824c71e0759a4f[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Tue Jun 8 15:30:56 2021 -0700

    Add function to make a disk harmonic mode basis
    
    The function accepts user-specified lists of m values and range of n values.
    Passed preliminary tests.

[33mcommit 0ad0ff87ca84d923ac6963fbadcb8c45717564df[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Mon Jun 7 16:09:29 2021 -0700

    Provide a draft version of function to alter DM1 using disk harmonic basis
    
    New code not yet tested.  Plan to first break up the function into more
    modular functions.

[33mcommit c0ed5c57b388d71e4f0920d695a3049c0b860d36[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Tue May 11 14:09:03 2021 -0700

    Create ipynb to demonstrate usage of hcipy disk harmonics

[33mcommit e6b8120705e048dc76290cf4f31f3b11601684ee[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Mon May 3 14:48:57 2021 -0700

    Save example for add a sleec feature
    
    See <diff_45f32_0bd2e.ext>.

[33mcommit 813cf5e40f3cd3186ba446471e1d006a11c3e97f[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Thu Apr 22 22:43:58 2021 -0700

    Update README.md

[33mcommit 98285c6a2ff4c8cf229092ce6ef111b895cbb28a[m
Merge: 5a10f10 2568331
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Jan 21 17:57:49 2021 -0800

    Merge branch 'master' of https://github.jpl.nasa.gov/pinchen/sleec

[33mcommit 5a10f10606ccbbeafd55288ce441675640515fca[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Jan 21 17:56:44 2021 -0800

    Make stylistic code corrections

[33mcommit 2568331aae2a9e147c9be313b94625b4a9360963[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Jan 21 17:56:44 2021 -0800

    Make stylistic corrections

[33mcommit 6e4b4b7dfdfd0f11200c9e5726ba6e7e0378353d[m
Merge: 1f38ad6 c104c17
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Thu Jan 21 17:44:20 2021 -0800

    Merge pull request #1 from dabdulah/master
    
    Add annular Lyot function

[33mcommit 1f38ad63a7b145bfb5c59ced98692fd8760fe774[m
Merge: 6490473 9de24a0
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Jan 19 17:59:44 2021 -0800

    Merge branch 'master' of https://github.jpl.nasa.gov/pinchen/sleec

[33mcommit 6490473b3cc75ccc13f221e69ca145ee1a9d596e[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Jan 19 17:57:07 2021 -0800

    Launch run <20210119t1756>, fix bug

[33mcommit 9de24a0bd3831463e3bf4b26d01ca7c1bb21fda0[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Tue Jan 19 16:59:19 2021 -0800

    Correct a typo

[33mcommit 763c2c91a30f1253c05d2596dd201b3de7d099d2[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Tue Jan 19 16:51:11 2021 -0800

    Remove obsolete ipynb files

[33mcommit 7971d9ad0837d7eb71ba61069a57cb580879c55e[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Tue Jan 19 16:45:27 2021 -0800

    Improve notes on Fourier transform of phase modulation

[33mcommit 059ae2b9578a7a39e666ff6af7001dbf84165ebd[m
Merge: f667ed0 edb1a06
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Tue Jan 19 14:18:44 2021 -0800

    Merge branch 'master' of https://github.jpl.nasa.gov/pinchen/sleec

[33mcommit f667ed036ba324bcecca3513ed7fd1acbc40dd14[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Wed Jan 13 23:40:03 2021 -0800

    Correct errors in mathematical equations
    
    Fourier transform of phase modulation.

[33mcommit edb1a06342a2169b267a9a76f428d12fb6382fe8[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Jan 13 18:05:20 2021 -0800

    Launch run <20210113t1802>

[33mcommit f8c2dc968d7a957e0f043649a54293a664b57414[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Jan 12 13:40:09 2021 -0800

    Launch run <20210112t1312>
    
    Also, fix minor bugs and update version to 0.2.0

[33mcommit fe2f7af4e95df8a4f38ef0db38b1c27f2deaf68d[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Jan 12 11:33:25 2021 -0800

    Fix typo

[33mcommit 6c0d06ed162b01774fc021255c10417ca154e457[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Mon Jan 11 18:44:38 2021 -0800

    Add plot and long-form options
    
    Can print header only and disable plotting now.

[33mcommit 8b2e69f0973d4f874d7cacaa79562acb8eedf0a5[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Mon Jan 11 18:30:39 2021 -0800

    Added planes shown in plotted ouput

[33mcommit 1621db50a908683ca0e2a69691eceef54a659e70[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Mon Jan 11 15:54:13 2021 -0800

    Update code slightly
    
    Remove extraneous `diameter` & `grid` arguments.  Updated parameter
    values in <sleec_pars.yml>.

[33mcommit bcbbcd7976c57ff778fa761491cf6db4fc0dc7af[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Mon Jan 11 14:51:49 2021 -0800

    Update README.md

[33mcommit 87be8b47a2cdfa663d9009e1004cfaa386a6ae5d[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Mon Jan 11 14:46:35 2021 -0800

    Enable usage of Slepian-Fourier basis for DM2 in EMCEE runs
    
    Ran <sleec_run> one iteration without runtime error.

[33mcommit a7ebdee50bf7ece4a106cc1f33358de47e4e74f7[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Mon Jan 11 14:05:16 2021 -0800

    Add function to update DM2 with opposing DM1 change
    
    New function not yet tested.

[33mcommit 257df2fb952269372e977011500f795c0e6ac9d7[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Fri Jan 8 18:41:33 2021 -0800

    Update README.md

[33mcommit 61ac9ae55dc932e113c31c7f4d21c73c14553773[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Fri Jan 8 18:40:56 2021 -0800

    Update README.md

[33mcommit 7a1e5a3b8aeb16aa2a2ccfe0c8a4141a33b6b47c[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Fri Jan 8 18:37:27 2021 -0800

    Add instructions for adding an optimization feature

[33mcommit 46f32b23334b3f5965ee6a05007bdc94df486e82[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Fri Jan 8 18:18:55 2021 -0800

    Enable usage of Slepian-Fourier basis for DM1 in EMCEE runs
    
    <sleec_run> runs without runtime error.

[33mcommit 0bd2ef219dc34f563a690a0cc15b2221bec67154[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Fri Jan 8 17:33:38 2021 -0800

    Rename `n2d` input parameter to `n2d_pupil`
    
    Making the parameter name more specific to accommodate using
    Slepian basis sets in planes other than the pupil.  Also
    removed extraneous parameter `n2d_coefs`.

[33mcommit ab88b6555492e50f57866a562ba7fab6f702b57d[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Thu Jan 7 13:09:11 2021 -0800

    Correct error in doc string re. `n2d`

[33mcommit c104c173ff2018fe13e26ed716376d4b513c8f84[m
Author: Daniel Abdulah <Daniel.Y.Abdulah@jpl.nasa.gov>
Date:   Wed Dec 30 17:09:45 2020 -0800

    add annular Lyot function

[33mcommit 5e4bc10f9ce9a6d94d2f04cff460e7bfcab90172[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Wed Dec 30 15:56:05 2020 -0800

    Update ipynb on 1-D phase modulation

[33mcommit 321ddf55e3f64bf10ca99bb2a263cb689a3562ec[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Tue Dec 29 22:03:30 2020 -0800

    Semi-quantitatively verify Slepian-Fourier DM1 update function

[33mcommit 1e5299a2c0ce3583370eec31d438cd00e41aacba[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Tue Dec 29 18:53:01 2020 -0800

    Preliminarily verify Slepian-Fourier DM1 update function

[33mcommit b6e0842f8d511288414c1a8e5423b97c7f3b3103[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Tue Dec 29 18:07:48 2020 -0800

    Add code to update DM1 usinig Slepian-Fourier basis

[33mcommit 4d946b4465657eaaf8cad2b89f526ad957bc0e83[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Wed Dec 23 22:39:04 2020 -0800

    Correct doc string for Slepian-Fourier generator function

[33mcommit 8d0590336cde71684597a1fc7d7cb81b06c921f7[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Wed Dec 23 21:42:51 2020 -0800

    Generate Slepian-Fourier basis set
    
    Output hcipy.ModeBasis object.  Qalitatively verified output using
    test routine.

[33mcommit 116ac2adf8ae5e682159f18f5f78458d3e823053[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Wed Dec 23 15:30:46 2020 -0800

    Update README.md
    
    Provide more detailed instructions on installing downloaded patches.

[33mcommit 240d128b4bfeeb943cedb0abb71b690163a19857[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Wed Dec 23 15:03:21 2020 -0800

    Update README.md

[33mcommit b03bad671fe61c30ef6d7e37012740b7cf52f806[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Wed Dec 23 15:02:08 2020 -0800

    Update README.md

[33mcommit 7b627efaefa945212ac3ba3217eb6849702bdc1d[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Wed Dec 23 15:01:11 2020 -0800

    Update README.md

[33mcommit e2977f0dd2de6d958f34250a7f341502987dc7bb[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Wed Dec 23 15:00:39 2020 -0800

    Update README.md

[33mcommit 75722365021999de9e802f7a8bfd71c6ffdf3854[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Wed Dec 23 14:51:31 2020 -0800

    Add instructions on installing patch fixes and update HCIPy prerequisite

[33mcommit 84c20d3d23ee4eacf9cbf3ac2ab67fbb0a760c56[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Wed Dec 23 12:52:37 2020 -0800

    Update html version of ipynb

[33mcommit 7f1e1dbe5b8dda25db19f7ff37753efd35cbca8b[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Wed Dec 23 12:49:02 2020 -0800

    Generalize ipynb ansatz to all real functions

[33mcommit b6e7b0f537280db47b2c05b118f433aa76e2a39c[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Wed Dec 23 11:52:52 2020 -0800

    IPYNB provides more accrately worded ansatz

[33mcommit 8f2b48c0d30655dbc912ef0b3de67d19f92f5e14[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Tue Dec 22 17:46:28 2020 -0800

    Ipynb provides mathematical ansatz on phase modulation

[33mcommit 2317b3099894b894c2b06d44a4222607cb5a6f04[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Dec 22 14:43:14 2020 -0800

    Git remove *.pyc in <sleec> node

[33mcommit ecf725af717aea54f22f0eb1946e2f0797db5f7f[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Tue Dec 22 14:36:54 2020 -0800

    Git will ignore .pyc files

[33mcommit ea1db0bb2402e87b51f3e490c04ac9c6d8fc7d4f[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Tue Dec 22 14:33:38 2020 -0800

    Git remove .pyc files

[33mcommit f4096202d73a3ab3fcd24ab21e5a9dfce62e2d30[m
Merge: 3097590 285118e
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Dec 22 14:24:36 2020 -0800

    Merge branch 'master' of https://github.jpl.nasa.gov/pinchen/sleec

[33mcommit 3097590ef2b5c60e4105e677446a4915211032ea[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Dec 22 14:24:12 2020 -0800

    Change which run instance to reconstruct

[33mcommit 285118e8b70249a82bd5e445be13fe56784a4fe1[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Mon Dec 14 17:24:14 2020 -0800

    Update README.md

[33mcommit bce475a5adddaece22c571fd1870f759f3222cde[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Wed Dec 9 17:39:48 2020 -0800

    In <misc> subdirectory, added .ipynb worksheets to simulate phase modulation.

[33mcommit 4dbd268efe858d85280c39951a6dcd3fba0d5357[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Sun Jul 26 15:36:42 2020 -0700

    <sleec_alter>: Found a bug.  <fpm_supergauss()> was missing a factor of
    f_num*wvl.  Fixed it.

[33mcommit 907ace253295d64a4778ba8e8eb10f9cb6ad8976[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Sun Jul 26 01:20:34 2020 -0700

    Launched run <20200726t0120>

[33mcommit 563868fb5da30550c58dac809e0e1a4d4943c324[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Sun Jul 26 01:15:54 2020 -0700

    Launched run <20200726t0114>

[33mcommit b400822fd3891e515b7cb91cf66b8c0acdd44fda[m
Merge: f75ea87 700c996
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Sun Jul 26 01:13:45 2020 -0700

    Merge branch 'master' of https://github.jpl.nasa.gov/pinchen/sleec

[33mcommit f75ea872cdc2a228c2d089b4b4bf961a4a7e98e4[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Sun Jul 26 01:12:39 2020 -0700

    Snapshot, no change since last commit

[33mcommit dc0123efd0f96d174656c70bdf0cfa0441d84bf1[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Sun Jul 26 01:09:50 2020 -0700

    Super-Gaussian FPM implemented.

[33mcommit 700c9962cba9a2a5978ea3f27cc9ab4327af067d[m
Merge: aa8b535 3a081dc
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Sat Jul 25 23:12:05 2020 -0700

    Merge branch 'master' of https://github.jpl.nasa.gov/pinchen/sleec

[33mcommit aa8b5357447668159e7a94e972b43adc25ec55d6[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Sat Jul 25 23:11:33 2020 -0700

    Snapshot

[33mcommit 6807f4dd95758eeac7d2161c92cd5e526e114c12[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Sat Jul 25 23:05:42 2020 -0700

    Launched run <20200725t2304>

[33mcommit 480443bb71e2fd5c45bc3060af5bb074824266d4[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Jul 23 23:07:49 2020 -0700

    Fixed minor bug.  Recent revisions in <sleec_run> and <sleec_pars.yaml>
    are working now.

[33mcommit d93d76f32206d771944ac13e7f6c891390149ed8[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Jul 23 23:04:48 2020 -0700

    Revised <sleec_run.run_emcee()>, untested

[33mcommit 07798a12718452a056b6743875fce044250c2e11[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Jul 23 22:39:26 2020 -0700

    Revised <sleec_run.ln_probability()>, untested.

[33mcommit e12fcb07a889120088a5695529dc96550a9012ea[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Jul 23 22:01:25 2020 -0700

    <sleec_run>: Revised <ln_prior()>, untested

[33mcommit 41be4267ea414d3e21926d560eb33dfa83d41c76[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Jul 23 21:47:16 2020 -0700

    <sleec_pars.yaml>: merged `bounds_xxx_pars` ('xxx' = 'pupil', 'fpm', or
    'lyot') into just `bounds_pars`.

[33mcommit e57c0ff52c8c11898155125c89f60ec9d30bce2f[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Jul 22 22:19:38 2020 -0700

    In <sleec_pars.yaml> and <sleec_run>:  Renamed variables `optmz_pupil_coefs`,
    `optmz_fpm_dia`, to `optmz_lyot_dia` to `bounds_pupil_pars`,
    `bounds_fpm_pars` and `bounds_lyot_pars`.
    
    Also, `bounds_pupil_pars` is now an Mx3 list, with M corresponding to
    the number of Slepian terms.  The 3 dimension corresponds to initial, min,
    and max values.
    
    Need to revise <ln_prior()>, <ln_probability>, and <run_emcee()> in
    <sleec_run> to be compatilbe with the new data structure.

[33mcommit 9fe22255b36b0c2e09c30f836f4845358f5307cc[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Jul 21 23:44:02 2020 -0700

    <sleec_alter>: Coded <fpm_supergauss()>, untested.

[33mcommit 6d51598bc5337023c3b5f8f83464119dc8e29776[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Jul 21 23:22:41 2020 -0700

    Reverted <sleec_pars.yml> to three-Slepian pupil, same as 2020-07-06 version.

[33mcommit 3a081dc0722c4aa96c0a4b2267bb6819a2d6c32e[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri Jul 17 22:01:26 2020 -0700

    Launched run <20200717t2200>

[33mcommit 15413b3fa4392630c0d24f85621baa1467d32e47[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Jul 16 22:09:48 2020 -0700

    Launched run <20200716t2200>

[33mcommit 87b02433d491e4402c7ed7582ea5f5766e4dd6ab[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Jul 9 18:13:52 2020 -0700

    <sleec_show>: Fixed bug in contrast plotting.

[33mcommit 3ae2a2cac654b35e9b0505370ec741f10d432a09[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Jul 8 00:15:33 2020 -0700

    <sleec_show>: Added <quad_plot()> and revised <reconstruct()> to improve
    plotting.

[33mcommit 55d17b2690735484ab956cddbd97ea683fdbee9d[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Jul 7 13:45:47 2020 -0700

    Updated <sleec_show.reconstruct()>.

[33mcommit 19db238d4134ea06187f7640d6596fd6f43fd8d4[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Mon Jul 6 23:02:39 2020 -0700

    Launched run <20200706t2300>

[33mcommit 56be7652b4c494ae7ba2e3618b197704431958af[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Mon Jul 6 16:42:18 2020 -0700

    <sleec_merit>: Catch when `m_0` == 0 and set it to eps to avoid `true_divide`
    runtime error. Launched run <20200706t1632>.

[33mcommit 789a9b827424dc8766723527c7e152ac0893f356[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri Jul 3 18:52:27 2020 -0700

    <sleec_run>: minor improvements.  Launched run <20200703t1851>.

[33mcommit 5cebb9a7082c3f7c7e55f8dc1340e2afe9ab8cb9[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri Jul 3 18:17:33 2020 -0700

    <sleec_run>: Revised to use new <sleec_merit.merit()>.  Run w/o error
    (not saving output).  Will launch recorded run <20200703t1818>.

[33mcommit c4ae419533998e4f5082a1e7da596c8db3a10764[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri Jul 3 17:28:33 2020 -0700

    <sleec_merit>: Removed `psf0_encircle` from argument list of functions; it
    is extraneous.

[33mcommit 017503f00eabfa3f1e686b9c1efaea110464332f[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri Jul 3 16:45:12 2020 -0700

    <sleec_merit>: Added `wavefronts` to the returns of <exposure_time()> and
    <merit()>

[33mcommit f814c57b0ae702c1a75695809cbadd1397ee1aeb[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri Jul 3 16:21:18 2020 -0700

    <sleec_merit>: Removed computation of `whole_contrast`.  It should save about
    0.17s of run time per merit call.

[33mcommit bee99fd70a27e40a587407a3dab0ad4efd3730f1[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri Jul 3 16:13:33 2020 -0700

    <sleec_merit>: Added new <merit()>.

[33mcommit d4bdab7b3c46549fa6cabc06c603670d0d10c8b2[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Jul 2 23:40:52 2020 -0700

    <sleec_merit>: Code cleanup.  <exposure_time()> checks out against .ipynb
    results

[33mcommit dd69350c1c3b17930c5cf8cac03414c17ac825d8[m
Merge: c20aeda f0b896e
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Jul 2 23:31:29 2020 -0700

    Merge branch 'master' of https://github.jpl.nasa.gov/pinchen/sleec

[33mcommit c20aeda83b9b87106acb98ba0e5421b79c74b000[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Jul 2 23:30:26 2020 -0700

    <sleec_merit>: Revised <exposure_time> to use encircled energy.

[33mcommit f0b896e7b44a84e6bab85a533c56f688ccc45180[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Thu Jul 2 23:29:16 2020 -0700

    Updated .ipynb parameter values for verifying <sleec_merit._test_exposure_time>

[33mcommit d34324647e2b8463178a81e51045b5d1c5c77a9b[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Jul 2 23:10:01 2020 -0700

    <sleec_merit>: revised <contrast()> to compute encircled energy (w/o FPM)
    instead of `penalty`

[33mcommit 2f8ddb852164743207b28aa7b1b1778dd66475ec[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Jul 2 22:47:55 2020 -0700

    <sleec_merit>: Replaced <psf0_variance> with <psf0_encircle>, which returns
    PSF encircled energy of a clear aperture.

[33mcommit f22de9da6b2ce0b780026ef426be5058d49760eb[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Jul 2 18:46:27 2020 -0700

    <sleec_merit>:  Editing <contrast()> in progress.

[33mcommit a277b05822cf4b4cfe8e4c57bd714d3d5123c2fb[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Jul 2 18:13:12 2020 -0700

    Verified <sleec_merit._test_exposure_time()> against
    <Merit_Function_Stark-etal.ipynb>.

[33mcommit 9d020f2012e6570554521f3483f6d5486cdd9fae[m
Merge: 9178b35 d154a2f
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Jul 2 18:12:28 2020 -0700

    Merge branch 'master' of https://github.jpl.nasa.gov/pinchen/sleec

[33mcommit d154a2feec1fb7e54da06dc34a1992e0d5637c3c[m
Merge: ed9a165 9f7931f
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Thu Jul 2 18:11:43 2020 -0700

    Merge branch 'master' of https://github.jpl.nasa.gov/pinchen/sleec

[33mcommit ed9a165ae823286f3497ca78d5617bab0cf0f1d4[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Thu Jul 2 18:10:19 2020 -0700

    <Merit_Function_Stark-etal.ipynb>: changed input parameter values for
    verifying <sleec_merit._test_exposure_time()>

[33mcommit 9178b35f13e95a2728e78a12d08528cb8d1430b6[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Jul 2 16:16:06 2020 -0700

    <sleec_merit>: Fixed bug in variance calculation

[33mcommit 7e033559115f0b6852ceb6dd0c65fc4833388edd[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Jul 2 16:03:27 2020 -0700

    <sleec_merit>:
    1. Renamed <psf0_footprint()> to <psf0_variance()>
    2. <contrast()>: Incorporated `psf0_variance` in `pentalty` calculation

[33mcommit 1dc7ec637b064414bf099366de774be25b396a64[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Jul 2 15:58:59 2020 -0700

    <sleec_optics>: Removed  from argument list of <Coronagraph>

[33mcommit 15a7487b968852526dd3454d3179f8510bdad017[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Jul 2 15:39:12 2020 -0700

    <sleec_merit>: Added <psf0_footprint()>

[33mcommit 415da828742831287ece28ccab2c56509a2aa956[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Jul 1 23:09:01 2020 -0700

    <sleec_merit>: Need to finish coding <exposure_time()>

[33mcommit d4a1ecd1897054e89cb5453f976abc7f9f9280b6[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Jul 1 22:43:28 2020 -0700

    <sleec_merit>: Eliminate option arguments `area_penalty` and
    `ouput_wavefronts` from <contrast()>. They are effectively always `True` now.

[33mcommit 24e41b1fcff81b2ddb2ad19643b7992ce226d970[m
Merge: 9f7931f db21a07
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Jul 1 22:18:07 2020 -0700

    Merge branch 'thruput'

[33mcommit 9f7931f8e9de1149f60813368c032a52f907ec3b[m
Merge: 0ba4ee4 6cab1b5
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Jul 1 22:17:04 2020 -0700

    Merge branch 'master' of https://github.jpl.nasa.gov/pinchen/sleec

[33mcommit db21a07460b024c0ae8bd95ab77fb89639ba6379[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Jul 1 17:56:54 2020 -0700

    <sleec_merit>: Added <transmit_all()> to compute a coronagraph's end-to-end
    transmittance/throughput.

[33mcommit 92e338bd942b7a956ba58ceecdb886cbfdf35a5d[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Jul 1 16:50:00 2020 -0700

    <sleec_merit>: Added functions to compute planetary, residual stellar,
    zodiacal, and exozodiacal count rates.

[33mcommit 6cab1b55b47d65dd885f087c9542b53b00dd17bd[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Wed Jul 1 16:35:28 2020 -0700

    Printed .ipynb to PDF

[33mcommit 94832c4f52e123872a52d7fa84f7949efa084c87[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Wed Jul 1 16:27:27 2020 -0700

    <Merit_Function_Stark-etal.ipynb>: Added computations of terms for baseline
    case.

[33mcommit 27b59fd108cfe09d494f7a19be4f185e612f1d01[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Tue Jun 30 14:09:55 2020 -0700

    Corrected typo in .ipynb file

[33mcommit 33cb6ab356dceed940ad1232eb3b483ca9f8a3e9[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Tue Jun 30 13:50:28 2020 -0700

    Added .ipynb write-up of Merit function based on Stark et al. (2014)

[33mcommit 63c60918c226417bb365deda97b6b969a86f23d2[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Thu Jun 25 14:18:36 2020 -0700

    Move laptop directories to new location.

[33mcommit 33f5451f99e72c26bfdaec592dc1e989929ad038[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Jun 16 23:14:40 2020 -0700

    <sleec_merit.contrast()>: Coded clean-up and doc strings.

[33mcommit d485296769878cf20b485fcc0664e32d25b7364f[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Jun 16 22:42:34 2020 -0700

    <sleec_merit>: Added code in <contrast()> to estimate footprint penalty in
    terms of second moment of intensity.

[33mcommit 57f5c493f4487ee99850525d13bb41bc4a80685e[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Mon Jun 15 23:44:06 2020 -0700

    <sleec_merit>: Working on revising `penalty` in <contrast()>.
    In <if area_penalty> block, need to properly define areas in terms of
    mean and variance.

[33mcommit e26eac9310658e476ed41e4ec24cb3351b77ceb9[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Mon Jun 15 22:19:25 2020 -0700

    <sleec_optics>: code clean-up

[33mcommit c911b7c940dd0672ef20613e5563fd5796e360fc[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Jun 10 14:37:50 2020 -0700

    <sleec_merit.py>: Added <fpm_transmit()>

[33mcommit a69613a17638cb6dbb55dea3d0e66c3563e64ab4[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Jun 9 18:25:18 2020 -0700

    <sleec_merit.py>: Added <lyot_transmit()>

[33mcommit 973cce692697e39b558073d8b56f3e813d0f7c4b[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Jun 9 18:06:54 2020 -0700

    <sleec_merit.py>: Added <pupil_transmit()>.

[33mcommit 11033b98cce828a4186035a46438148c44c3343b[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Jun 9 17:21:16 2020 -0700

    Replaced variable names `pupil_diamtr` and `dm_diamtr` by `pupil_diameter`
    and `dm_diameter` in all code.

[33mcommit a0b02045d563d0be6b33c46f1e9c69242d3639df[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Jun 9 17:12:26 2020 -0700

    snapshot

[33mcommit 7d1c66630f4a72d7e200f18a1443d7e73e37e7a2[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed May 27 16:32:39 2020 -0700

    <sleec_merit.py>  Renamed <throughput_annul()> to <throughput_annls()>

[33mcommit fcb89a3d3d2008f7a779c8e1db9f7feabc11f08a[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Mon May 18 17:58:56 2020 -0700

    <sleec_merit.py>  Renamed <throughput()> to <throughput_rand()>.  Brought
    back old <throughput()> (commit 501ee87) as <throughput_annul()>, with
    modified return options.

[33mcommit 4ab9aca71282962773beb60a68e5739264e92a1b[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Sun May 17 17:42:18 2020 -0700

    Launched run <20200517t1741>

[33mcommit daf7132f8b84b4e1f13bcb268237f623e331b43d[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Sun May 17 17:37:20 2020 -0700

    Snapshot.  Run <20200516t0007> was launched in this state.

[33mcommit 4d33da9ec2352371bbe99a99bc5de59ab396aa60[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Sat May 16 00:00:49 2020 -0700

    Revised <sleec_run.py> and <sleec_show.py>.  Ran <sleec_run.py> a few
    emcee iterations w/o runtime error.

[33mcommit 764c3d0be8765896d95c0cf56deae86209c5a91e[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri May 15 23:36:41 2020 -0700

    <sleec_merit.py> New <throughput()> computation by generating a planet
    randomly positioned in the dark zone (instead of illuminating the entire
    dark zone) is working.
    
    Need to revise other modules, <sleec_run.py> and <sleec_show.py>, to make
    them compatible with the new interface.

[33mcommit 0ba4ee42f0abacd8457ee4ad51da8752db2ce3ae[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri May 15 21:30:11 2020 -0700

    Snapshot

[33mcommit 25dd81040856be1ed57737a5e78bfa5ac53c40e8[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri May 15 16:59:09 2020 -0700

    Added option to apply weights in contrast calculation

[33mcommit 4a1bc9d73e8226b57f43731db145be9c525598c6[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed May 13 17:29:34 2020 -0700

    Minor tweaks

[33mcommit 828c857e70b972f460394ced1932010f4a6414b2[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed May 13 17:01:18 2020 -0700

    Snapshot

[33mcommit 501ee87c1e1dfc97c3f34b104c98e5f99611dd32[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed May 13 16:37:20 2020 -0700

    <sleec_input.py>: Initializing coronagraph with Slepian pupil is working
    (tested) now.
    Using `oversample` = 12 instead of 6 reduced the four-fold artifact in the
    final (contrast) image by a lot.

[33mcommit 3005bd0d57e9af4249d974ffafc1bcc83e0b2041[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed May 13 16:19:22 2020 -0700

    Edited <sleec_show.py> to work w/ new return scheme, no testing done yet

[33mcommit 505eb91a5193b1594b03309b6ed7e096577771e3[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed May 13 16:18:09 2020 -0700

    Edited <sleec_run.py> to work w/ new return scheme, no testing done yet

[33mcommit 17c0d6432c2ce920af4ec327c1a62ffc1f89fa4f[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed May 13 16:15:45 2020 -0700

    Edited <sleec_merit.py> to work w/ new return scheme, no testing done yet

[33mcommit bcb5db777f58e37f472ed503be1e839ce06c6f29[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed May 13 14:47:04 2020 -0700

    <sleec_input.py>: Revised & tested <_demo()> to work w/ new return scheme.

[33mcommit de7d87e2209bec87d93be4999b1624df5acc98f4[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed May 13 14:36:49 2020 -0700

    <sleec_input.py>: <import_pars()> now returns both the coronagraph and pars.
    Revised and tested <sleec_alter.py> to work this new return scheme.

[33mcommit e245cf63428f788981ead9f8e8feb6d26cedc002[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Mon May 11 17:34:09 2020 -0700

    Revised code to add <slepian_n2d> option for initializing coronagraph pupil.

[33mcommit 91ea787faf407d57a20c0ff44c01a4b06a259e8e[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Mon May 11 12:43:38 2020 -0700

    Snapshot

[33mcommit edf67855a4f078ceb86296747682820c4deb9c6d[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri May 8 18:00:30 2020 -0700

    Launched run <20200508t1759>.  Same as previous (20200430...) run,
    except inclusion of `grid_oversize` and setting 5000 iterations.

[33mcommit 4276052378b3be35434297b6b1c24e4aa112d8b6[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri May 8 17:38:21 2020 -0700

    <sleec_input.py> and <sleec_pars.yml>
    Added 'grid_oversize' configuration parameter to provide margin between
    the square boundaries of the pupil grid and the pupil aperture's
    circumference.  This removed the four-lobe intensity pattern in the final
    image.  One remaining possible artifact is an azimuthal wave around the
    10th airy null ring in the final image.  Applying the margin in the focal
    plane did not remove this artifact.

[33mcommit 36be8f76a03a5c71792232f67154c7f39ea8151e[m
Merge: a917a1e dd18d02
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed May 6 17:47:29 2020 -0700

    Merge branch 'master' of https://github.jpl.nasa.gov/pinchen/sleec

[33mcommit a917a1eaa37074728e634342cc7887df16ba8d4e[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed May 6 17:47:17 2020 -0700

    Updated <sleec_show.py>

[33mcommit dd18d0247e7bede90352108e1b72be3fda6d0feb[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Wed May 6 17:37:48 2020 -0700

    Added <super_Gaussian_demo.ipynb>

[33mcommit 593eabed91f647f584a008d57d3fba2ffed5875e[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Apr 30 00:45:00 2020 -0700

    <sleec_merit> Had to add several patches to fix the chain of true_divide
    problems.  Launched run <20200430t0032>.  37 iterations so far; hopefully
    the problem is fixed.

[33mcommit b40ef8f0ed763acf5e8c90a3d2365ea8f62db547[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Apr 29 23:12:26 2020 -0700

    <sleec_merit> Previous patch did not work.  Added patch of setting
    `psf_core_area` = 1 when `psf0` is invalid, and set `psf0` = computer
    eps (i.e. reverted previous patch).
    Launched run <20200429t2311>.

[33mcommit ad63d1f0bf1b44125df46307cef1b62b28586823[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Apr 29 22:52:01 2020 -0700

    <sleec_merit> Run crashed due to true division error caused by `psf0` again.
    Patched it by setting `psf0` to ten times computer eps.
    Launched run <20200429t2251>

[33mcommit 379aa6bcf61bda235081c46b7183b2799475d551[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Apr 29 22:46:29 2020 -0700

    Updated <README.md>

[33mcommit b8dc01d6848dc08ac4b432cef07f07a9a4435b19[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Apr 29 22:40:55 2020 -0700

    <sleec_alter> Added circular aperture in <pupil_slepians_r()>. This
    minimizes artifacts of square PSF and speckle field w/ 4-fold symmetry.
    Launched run <20200429t2239>.

[33mcommit 25d662c46de7e1dc7d38c7bf60a07e3c7b7a4a1e[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Apr 29 22:19:30 2020 -0700

    <sleec_show.py> Made minor improvements.

[33mcommit b5fdd8493df81b2b22a019cb67cffccaae101cdb[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Apr 29 18:33:02 2020 -0700

    <sleec_alter> Fixed serious error:  <fpm_dot()> did not actually
    update the coronagraph.
    <sleec_alter> and <sleec_input> Fixed another error:  Field stop diameter was
    a factor of two too small.
    
    Launched run <20200429t1836>

[33mcommit 95b650e197e29dc365eb915b46881cd25bda04ac[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Apr 29 17:49:18 2020 -0700

    <sleec_run.py> Reinstated `non_negative` as an argument for
    `ln_probaility()`, to make <sleec_show.reconstruct()> work.
    Also, <sleec_show.reconstruct()> now plots images with
    pupil diameter units in pupil planes and airy units in image planes.

[33mcommit cb7a7dfa4d99cab20a248487cd238d9679b5bde2[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Apr 29 15:51:18 2020 -0700

    Moved <non_negative> option to become a prescribed parameter in
    <sleec_pars.yml>

[33mcommit e0a41b3950de99c74867751f454e988c4c9ff683[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Apr 29 15:29:41 2020 -0700

    <sleec_show.py> <reconstruct()> now works.  Will refine plotting and
    functionalities.

[33mcommit bbea8a55fee208d82cfe90ad66d9ce7906b3cbcd[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Apr 29 00:05:03 2020 -0700

    EMCEE is running w/o runtime error.  Run prefix = 20200429t0004.

[33mcommit 962145138ff63def29ffa26cde3633ae2155bd59[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Apr 28 23:55:53 2020 -0700

    Re-instated code for  option.

[33mcommit 90c4c0af31b7e273cb34e10114dd0438b63284cc[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Apr 28 23:39:18 2020 -0700

    Need to re-instate  option

[33mcommit c171ea923ef07bf3893ff899b3495ef681679c8a[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Apr 28 23:14:33 2020 -0700

    <sleec_show.py> Minor tweaks

[33mcommit cb5740be41d6deade05a6bf9072f62585e72ae46[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Apr 28 22:33:38 2020 -0700

    Fixed typos in <sleec_run.py>

[33mcommit f5db06ae4418fc5369ae11b39da74f3adf39081c[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Apr 28 22:22:35 2020 -0700

    <sleec_run.py>  Inserted code to copy <sleec_pars.yml> and <slepians.pkl>.
    Need to test code.

[33mcommit 30f49aa9f1f114a302e496d20a0117ac6f5368f1[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Apr 28 15:56:21 2020 -0700

    Coded option to reconstruct coronagraph and report outputs for a given sample
    from an EMCEE chain.  Need to revise code to save <slepian.pkl> and
    <sleec_pars.yaml> when a run starts.

[33mcommit 82f11a55892ed2b61f9f9450eb38ad5b3c9361ef[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Apr 28 14:51:30 2020 -0700

    Purged  option

[33mcommit 0184e69fac71a3f165255d4209f068dedf424f50[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Apr 28 14:04:04 2020 -0700

    <sleec_merit.py>  Fixed a bug in previous patch.  Launched run
    <20200428t1402.hdf>.

[33mcommit 67a972f1e9a0b06194ff04c39e1f5fdb329ee0ee[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Apr 28 13:50:49 2020 -0700

    <sleec_merit.py> EMCEE runs sometimes crash due to invalid division error
    associated with `psf.max()` in <sleec_merit.py>.  Patched it by setting the
    value to `np.finfo(float).eps` if `psf0` = `psf.max()` is not finite or is
    zero.
    
    Started run <20200428t1347.hdf>.

[33mcommit 7b689c0760228e03096457eab10068f51c4736f5[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Mon Apr 27 23:47:18 2020 -0700

    <sleec_run.py>  Re-instated prior limits on Slepian coefficients, hoping
    this would avoid divide-by-zero runtime crash.  Also, revised
    `ln_probability` function to avoid merit computations if prior rejection
    occurs.

[33mcommit c7b4c673f1a97e2ee3e5bc97c48bb908812c1d14[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Mon Apr 27 21:31:56 2020 -0700

    <sleec_show.py>  Added printing of throughput, contrast, planet_core_area,
    and area_penalty values for the top five samples of the EMCEE chain.

[33mcommit b50b029fe50ca66b91e3f867698472cc88b1981f[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Mon Apr 27 18:16:03 2020 -0700

    <sleec_run.py>  Fixed bug in limitting `pupil_field` values and added
    option to replace negative values with zero.

[33mcommit 7dc0ff91513ebe963a87743d398aa26485110d30[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Mon Apr 27 11:43:12 2020 -0700

    Added emcee 'blobs' functionality to <sleec_run.py>.
    Saving lists to HDF5 causes an error, so only saving scalars right now, which
    excludes saving `throughputs` and `contrasts`.
    In the future, saving list elements individually should work.

[33mcommit e4a152450d706e13237f8fc71a428b485d391802[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Mon Apr 27 01:41:16 2020 -0700

    Trial emcee run with 1000 iterations completed.  Corner plot shows all 5
    pupil Slepian coefficients are highly correlated.  Will launch a new run with
    only 1 coefficient.
    
    Added <sleec_show.py>.

[33mcommit 03a5ae409b1cdb171bdc960fa5654cdf9bf3e3d4[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Sun Apr 26 00:40:45 2020 -0700

    Updated version # to 0.1.0 (from 0.0.0)

[33mcommit f10d26206b8500080b4ec3972a7f8bf2dd010aa8[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Sun Apr 26 00:34:54 2020 -0700

    Updated README.md

[33mcommit ba9a91e2f1a7ca7aeed04059d58bd5d75d00bd1b[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Sun Apr 26 00:30:41 2020 -0700

    Multiprocessing works now.  Thanks to E. Por's hcipy patches for pickling issue.

[33mcommit acc80c00be2e29414db48ee49282afeb11f659f7[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Sat Apr 25 20:15:44 2020 -0700

    Minor tweaks

[33mcommit 117b7dc353274a6b6df10d79c1bf5c6344df95f2[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Sat Apr 25 19:58:11 2020 -0700

    <sleec_run.py> Minor tweaks

[33mcommit 45c0632d934d8307ce22bfd53985fa43e2350be4[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Sat Apr 25 19:19:50 2020 -0700

    New approach of using a pupil basis set of Slepians with m=0, k=1, and a
    range of N2D values is coded.  Running w/o runtime error.  Still need to
    verify and validate.

[33mcommit 9c5f03e6bf22bc2068beda3c7ed998f5caa8f17e[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Sat Apr 25 18:57:45 2020 -0700

    <gen_slepians.py> Revised for N2D ranging

[33mcommit c005240efc79419b6698f09547dc5deffdb69af9[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Sat Apr 25 17:57:41 2020 -0700

    <sleec_alter.py> Changed `pupil_slepian_r()` to accept range in N2D only

[33mcommit 306d22583a2f845e52d735906df35c382bb7224d[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Sat Apr 25 17:01:13 2020 -0700

    <sleec_run.py> Added  option to

[33mcommit fba20b08a491cb7594b76c90c35db5de271fcd40[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Sat Apr 25 02:23:10 2020 -0700

    Testing run_emcee.  Needs improvements

[33mcommit dcce541abac35bab687059a5064fa06813ba8142[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri Apr 24 19:32:17 2020 -0700

    Updated README.md

[33mcommit 7c3bd230987ebbc8e6b32ad4569e942928305f7b[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri Apr 24 19:28:22 2020 -0700

    <sleec_run>:  <run_emcee()> now runs w/o runtime error.

[33mcommit 1f33d853eb000685a6c7d2613fe3e8224c824fce[m
Merge: cec14fc 0d3435c
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri Apr 24 00:46:05 2020 -0700

    Merge branch 'master' of https://github.jpl.nasa.gov/pinchen/sleec

[33mcommit cec14fc24a7464fd0ce857de87718a2ddcf20267[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri Apr 24 00:45:38 2020 -0700

    <sleec_run.py> coding in progress

[33mcommit 0d3435c24bd16545c5ce3da5143b9b9f6e4dd9e4[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Thu Apr 23 00:44:02 2020 -0700

    Update README.md

[33mcommit 25e5ac881591ade133493e6592e5183fb525acb9[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Apr 23 00:40:47 2020 -0700

    <sleec_alter.py>:  Added <pupil_slepian_r()> function.

[33mcommit a4f8578cb8cfeacf845ab83feac8622df373672d[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Apr 22 18:09:17 2020 -0700

    Updating <sleec_alter.py>

[33mcommit 498a29bb0f746e5b665d7af60d90af5127a33538[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Apr 21 23:10:16 2020 -0700

    <sleec_merit.py>:  Added <merit()> function.

[33mcommit aeaf08a04f666d82ed64a50a34f45e41102fdfcc[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Apr 21 21:55:48 2020 -0700

    <sleec_merit.py>: Added <contrast()> function.

[33mcommit a029695273f06fe4b4e7771dbd7c7d3f1597b4df[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Apr 21 16:44:36 2020 -0700

    <sleec_merit.py>:  Added option to insert planet and output number of planet
    core pixels in <throughput()>.  Overall throughput with and without planet
    agree within 0.4%.

[33mcommit 32eb5ac54a47ce88a426fbf78774bb6e432fc2e7[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Apr 21 00:29:16 2020 -0700

    <sleec_merit.py> in the process of coding the cost function.

[33mcommit 383247d639118ccaac5b612c0b68d098381b7ad0[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Mon Apr 20 22:57:46 2020 -0700

    Eliminated <dark_field> attribute, which is no longer used, from
    `sleec_optics.Coronagraph` class.

[33mcommit 1bb973cfcc37468bf5dce7ed19cd6722d96f9703[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Mon Apr 20 18:58:42 2020 -0700

    <sleec_merit.py>. Added <wavfrnt_off_axis()> function.

[33mcommit a3274067ea3823b1c51194f10a1c3d6c9111835d[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Mon Apr 20 10:11:06 2020 -0700

    <sleec_merit.py>:  Initial version of `throughput()` is working.

[33mcommit 27cdbf6a20bf5db76640c325a3858f60f75bec17[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Sun Apr 19 18:17:35 2020 -0700

    <sleec_input.py>:  Replaced <hcipy.make_focal_grid_from_pupil_grid> with
    <hcipy.make_focal_grid>.  This simply removed depracated <hcipy> usage.

[33mcommit 62331a8e2952680a26bfb0fd72e9dc7548d607cb[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Sun Apr 19 16:21:28 2020 -0700

    <sleec_merit.py>:  `wvfrnt_annls_cnj()` now working.

[33mcommit d7d7adc37448aa631efc256741df1767886b23d6[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri Apr 17 18:36:39 2020 -0700

    <sleec_input.py>:  Added doc strings.

[33mcommit 38cc1e814b04877b740a7abbe2e4b0f55cb1af10[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Thu Apr 16 19:11:49 2020 -0700

    Update README.md

[33mcommit 5fafb35c96e5b9d6fd229afa4a86ece71f81a94e[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Apr 16 19:00:45 2020 -0700

    Reorganized code into more modules:
    <sleec_optics>:  Remains the same
    <sleec_input>:  Now has just a function to to import parameters from the
        <sleec_pars.yml> file and initialize a coronagraph object.  `_demo()`
        demonstrates the usage.  It now runs w/ the current yml file; have not
        tested for other cases.
    <sleec_optimize>:  Functions for computing merit/cost values were moved out
        wof <sleec_input> into here.  This module is in a unorganized state.

[33mcommit 6ce911e1f2e7f8053e78b3ed45fa5794069b71d4[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Mon Apr 13 14:01:36 2020 -0700

    Edited <README.md>

[33mcommit c3515c3374774e066b97ea94e94355d90774e5b4[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Apr 9 17:09:21 2020 -0700

    Added <sleec_pars.yml> (modeled after J. Trauger's <sleecParms.yml>) as the
    user interface for specifying initial parameters.
    Plan forward:
    1. Integrate update functions from <sleec_input.py> into <sleec_optics.py>
    `Coronagraph` class.
    2. Write a method to parse the yml file in <sleec_input.py> and clean up
    the rest of <sleec_input.py>.
    3. In <sleec_optics.py>, write methods to compute
        a. single annular contrast
        b. mean contrast in dark zone
        c. thruput
    4. In <sleec_input.py>, write cost function
    5. In <sleec_input.py>, write function to do <emcee> optimization.

[33mcommit 2de6412e261e0bce0c2471acaec758a666960197[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Wed Apr 8 18:40:20 2020 -0700

    Update README.md

[33mcommit 0ab4a556da6c6bff0e5039de3eabc0f342e75704[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Wed Apr 8 18:38:04 2020 -0700

    Update README.md

[33mcommit 04ef4085bbbae4cb3d98a4935f59e76d736f0dc6[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Wed Apr 8 18:31:46 2020 -0700

    Update README.md

[33mcommit bd793a582d0c0974b722c87264a2c5760c9c84b0[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Apr 8 18:24:34 2020 -0700

    Imported turbo color map

[33mcommit dc7941b9a7788c8a99a635d676c3c32b8b0d415f[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri Apr 3 22:25:10 2020 -0700

    <update_dm1_surfc()> and <update_dm2_surfc()> functions in <sleec_input.py>
    working now.

[33mcommit fbbfa584a0e78dfe27757fe28d6cc6c4ed722da3[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Apr 2 23:25:13 2020 -0700

    Removed .swp file

[33mcommit 0cef1765613a7f9259bb21a4dba2b2c5fffc1530[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Apr 2 23:21:17 2020 -0700

    Removed obsolete <hcipy> and <hcipy_local> directories.
    Added exception trapping to fix bug when `dm_surfc` is `None`.

[33mcommit ba8d2928acf8e02001a6a36a8360c4b731fe311c[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Thu Apr 2 22:37:27 2020 -0700

    Testing running  <sleec> on a new machine (JPL laptop w/ Ubuntu in vbox).
    Pip installed latest <hcipy> release v0.3.1.
    Kept getting <imsave> error that occurred with old <hcipy> when running
    <sleec_input.py> or <sleec_optics.py>.
    
    Diagnosis:
    ---------
    In director <~/sleec_out/> the .pkl influence-functon file has the old hcipy
    <hcipy_local> imbedded in it.  So in <sleec_optics.py>, when <load_element()>
    loads the .pkl file, <hcipy_local> is imported/loaded.  This raises the
    <imsave> error.  Although when running <sleec_input.py> as __main__, the
    <demo_dmsurfc> function does not use influence functions, sheer loading
    of <sleec_input> initiates a chain of events leading to the error.  The key
    step is the <ini_corngrph()> ftn, its default argument value of
    <default_args()> leads to loading the .pkl file.  Simply getting rid of
    the default value fixed the bug.
    
    Now sleec should run on a new machine w/ hcipy 0.3.1 or later installed.

[33mcommit 7ad2b7e4844b0e0e63a39fe95ba276bea56f96c6[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Apr 1 23:28:56 2020 -0700

    Added option to directly input DM surface figure, instead of using only
    vmap.  New code runs w/o runtime error and output plots of intensities
    in optical planes look OK, with flat DMs.  Need to test w/ DM shapes and add
    update-DM-surface function/option.

[33mcommit 1ca62cdac1537efb318ae4d0efe100974beb469e[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Apr 1 00:01:37 2020 -0700

    <slepian_simons.py>:  Corrected some doc-string errors.

[33mcommit 4d768f0fb7e33a463b275d04f249910e37b8c5e6[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Mar 31 18:40:21 2020 -0700

    Wrote script to compute Slepain functions.  See <slepian_simons> subdirectory.

[33mcommit 4914529da6bf26cacff7516de5667387a5bc3974[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Mon Mar 23 17:14:26 2020 -0700

    Added <update_dm1_vmap()> and <update_dm2_vmap()>.  <_test_dmX_vmap()> show
    expected results.

[33mcommit bc6bc0dc1d291e7588afb91032da1ecc654b033e[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri Mar 20 15:00:50 2020 -0700

    <README.md>: Re-worded to-do list.

[33mcommit 66a6f2400d45b254f0ad422b75a818b3953e738e[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Fri Mar 20 14:55:35 2020 -0700

    Updated to-do list in README per teleconference.

[33mcommit d0d2a287ac1013398a6d7941308ad7bd3c2269f7[m
Merge: 4f411f7 10cc1da
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Feb 20 21:17:21 2020 -0800

    Merge branch 'master' of https://github.jpl.nasa.gov/pinchen/sleec

[33mcommit 4f411f71dd533d922bd11662f07c22bd61e0a038[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Feb 20 21:09:12 2020 -0800

    <sleec_input.py>:  Debugged <mean_contrast()>. The cause of error had
    nothing to do with "field-stop artifacts."  The revised version yields
    fairly constant average contrast:  2.10E-7 @ num_airy=256, 2.21E-7 @
    num_airy=64, 2.67E-7 @ num_airy=32.

[33mcommit 10cc1daf7cf2ed50a3a08c881cea832f37688ceb[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Feb 20 21:09:12 2020 -0800

    <sleec_input.py>:  Debugged <mean_contrast()>. The revised version yields
    fairly constant average contrast:  2.10E-7 @ num_airy=256, 2.21E-7 @
    num_airy=64, 2.67E-7 @ num_airy=32.

[33mcommit abc93b7d6744e915c72e2c5e94bb5ea819772cb1[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Feb 19 22:30:06 2020 -0800

    <sleec_input.py>:  Debugging <mean_contrast()> result decreasing w/
    increasing <num_airy>.  Revised <mean_contrast()> so I can plot intermediate
    results.  The dark field size and peak intensity remain constant w.r.t.
    <num_airy>.  However, faint airy rings become more intense as one
    decreases <num_airy>.  I suspect that the field stop might be causing
    diffractive artifacts.

[33mcommit 0b7e1978f3b264bc2ad44457978cca95449f1140[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Feb 4 17:05:36 2020 -0800

    Edited <README.md>

[33mcommit fc317a6f6737aca8bf95a23bad09e0e7bd418da2[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Feb 4 16:53:57 2020 -0800

    Edited 'README.md'

[33mcommit 13e9d246375589ff88e94830147e8fdd79ac8b62[m
Merge: 6875892 b2b436e
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Feb 4 16:40:42 2020 -0800

    Merge branch 'hcipy'.  Resolved merge conflicts in <sleec_optics.py> and
    <sleec_input.py>.

[33mcommit b2b436e4d98aa7e3bd78acf3a5100f299eb3225d[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri Jan 24 09:59:59 2020 -0800

    <sleec_input.py> working w/ newest HCIPy ver now.

[33mcommit 1a10adac56b24a0a42d5ac01abfaec4155d05173[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Jan 23 17:48:29 2020 -0800

    No code change from last commit.  Emiel Por fixed the bug in HCIPy today.
    After installing the fixed release, <sleec_optics.py> works now.

[33mcommit 68758921f47199602e80496c7fe1cda05f3be3bb[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Jan 23 17:42:23 2020 -0800

    Commit for switching to hcipy branch

[33mcommit 3a6bdae20828c2f2a38c206a6b8c7a3fc1a54026[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Jan 22 21:59:09 2020 -0800

    <sleec_input.py>:  Optimized w/ <oversmpl> = 4 instead of 8.
    Final avg contrast = 1.2e-8. f value = 2.5e-8.
        FPM dia. = 5.30 wvl/D, Lyot dia. = 0.706D.

[33mcommit 1b4a5c517496dbf15fa4714456f7445c8b60f76b[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Jan 22 20:50:26 2020 -0800

    <sleec_input.py>:  Optimized w/ <cost_ftn>, which basically minimizes
    integration time.  Final avg contrast = 1.2e-8. f value = 2.4e-8.
    FPM dia. = 5.37 wvl/D, Lyot dia. = 0.702D.

[33mcommit 1dfaad81226bcebb7c47f13818cb6dc1380e7703[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Jan 22 20:40:23 2020 -0800

    <sleec_input.py>:  with <num_airy> = 2*48 and optimizing both FPM and
    Lyot-stop diameters simultaneously for avg. contrast, got 4e-9 avg. contrast.
    FPM dia. = 5.86 wvl/D, Lyot dia. = 0.574D
    See <~/sleec_out/cost_ftn3.log>.

[33mcommit f3dd6a2554647f5664d9f0343586fd12fbfe1fde[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Jan 22 20:03:22 2020 -0800

    <sleec_input.py>:  Added <cost_ftn_x()> and <_try_mystic_xxx()> functions for
    optimizing one parameter at a time to better understand what the optimiation
    is really doing.
    Found that using too small a value for <num_airy> caused artifacts (a dark
    ring in camera image at small angular radius, depending on the FPM and/or
    Lyot-stop diameter).  Changed <num_airy> from 28 to 54, and the optimzed Lyot-
    stop diameter changed from .36 to .73w .

[33mcommit 478473c466a7e0f66c12de002e450f143ea2bd78[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Jan 22 15:27:08 2020 -0800

    After getting rid of <wavelength_0> variable in a few places, the code runs
    with HCIPY v0.3.0.  But I now get the following error (excerpt):
    --EXCERPT OF ERROR MESSAGE--
      File "/home/pinchen/sleec/source/sleec_optics.py", line 159, in forward
        return self.prop.backward(wavefront)
      File "/home/pinchen/anaconda3/lib/python3.7/site-packages/hcipy/optics/optical_element.py", line 761, in res
        instance_data = self.get_instance_data(None, wavefront.grid, wavefront.wavelength)
      File "/home/pinchen/anaconda3/lib/python3.7/site-packages/hcipy/optics/optical_element.py", line 633, in get_instance_data
        self.make_instance(instance_data, input_grid, output_grid, wavelength)
      File "/home/pinchen/anaconda3/lib/python3.7/site-packages/hcipy/propagation/fraunhofer.py", line 38, in make_instance
        instance_data.fourier_transform = make_fourier_transform(input_grid, instance_data.uv_grid)
      File "/home/pinchen/anaconda3/lib/python3.7/site-packages/hcipy/fourier/fourier_transform.py", line 162, in make_fourier_transform
        if input_grid.is_separated and input_grid.is_('cartesian') and output_grid.is_separated and output_grid.is_('cartesian') and input_grid.ndim in [1,2]:
    AttributeError: 'function' object has no attribute 'is_separated'
    -----
    E-mailed Emiel Por for advice.

[33mcommit 6342c5f19613e1c4ed33fa34356e06cee7bbdf96[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Jan 22 14:04:29 2020 -0800

    E. Por, author of HCIPY, informed me that v0.3 just came out yesterday, and
    it should not have the backward-propagation issue, which raised "Output grid
    unknown, propagate forward first."  I updated to HCIPY v0.3 and copied the
    original <sleec_optics.py> and <sleec_input.py> files from the <master>
    branch.  Will now try how they run in HCIPY v0.3 environment.

[33mcommit d132391c37e77078ddd267d0098faaa553aceed3[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Jan 22 14:02:02 2020 -0800

    Snapshot for switching to 'hcipy' branch

[33mcommit 355997062517c3ad6b3ac7eb6f4c82f66c2948b4[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Jan 21 17:31:20 2020 -0800

    Branch <hcipy>:  Abandoning trying to make our sleec code compatible w/
    hcipy v0.2 (latest release).  The added complexity of the new release might
    not be what we want.

[33mcommit 4c5eac20355b348fe2158be695555900c062671d[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Jan 21 13:47:16 2020 -0800

    Ticket closed.  Verified in <sleec_input._trymystic()> that <coro.fpm_field>
    contains the FPM amplitude transmission data.

[33mcommit a9be145260a6ce683edde745529f3523a57a053f[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Jan 21 11:53:34 2020 -0800

    Committing so I can switch back to master branch

[33mcommit bd6a6d100bec4c2acbdd8ad2f6fa2342d3a575b8[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Mon Jan 20 22:28:22 2020 -0800

    <sleec_optics.py>:  Working on making our code compatible with <hcipy> v0.2.
    Need to eiher update <FresnelBack> class to avoid error stating "Grid
    unknown.  Propagate forward first."  Or figure out how to use HCIPY's native
    <FresnelPropagator> class.
    Added <verbose> attribute for <Coronagraph> class in <sleec_optics.py> and
    <ini_corngrph> function in <sleec_input.py>.

[33mcommit a13b81db8ba8f6edc81b6a0186e17a33c29e80fc[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Sat Dec 28 23:47:54 2019 -0800

    <sleec_input.py> runs on <sleec> node, no code change from previous commit.

[33mcommit fd6dc1f432457e724885544c60806673bb4f01a7[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Sat Dec 28 23:08:48 2019 -0800

    <sleec_input.py>:
    The version in the previous commit was erroneous; it concatenated both of
    the conflicting versions of <sleec_input.py>.  I deleted the older (upper)
    half of the file.  Also, in <try_mystic.py>, I reverted to using
    <cost_ftn()> instead of <cost_ftn_2()>.
    
    Note:
    The current <hcipy_local> package is incompatible with <scipy> 1.2.0 or
    higher.  The problem is a occurrence of <from scipy.misc import imsave>
    statement.  <scipy> deprecated <imsave> and removed it from ver. 1.2.0.

[33mcommit 562ec40ffe2d59499fe21b0c79af9dd2393133fe[m
Merge: 90bafa4 f9c7a48
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Sat Dec 28 21:25:15 2019 -0800

    Edited <sleec_input.py> to resolve pull conflict

[33mcommit 90bafa4e9abb979368e3e121816f58e69eb4f8f4[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Sat Dec 28 21:12:13 2019 -0800

    Snapshot before pull

[33mcommit f9c7a481157f2e767a3bbccbd70886bb2d0895b1[m
Merge: 951ae20 6b8e37f
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri Dec 13 11:06:35 2019 -0800

    Merge branch 'master' of https://github.jpl.nasa.gov/pinchen/sleec

[33mcommit 6b8e37fb02d0aeb565536ef2c93d34f1b4d2cc79[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Fri Jul 12 17:32:56 2019 -0700

    <sleec_input.py>:  <_try_mysitc> optimization with field_stop
    working.  Varying initial values (of fpm width and Lyot-stop
    diameter) changes solution (i.e. degenerate solut'n space).

[33mcommit 951ae20a090c4bb408b97e41febde77ced95f5f7[m
Merge: 342e3d3 12d8dcd
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri Jul 12 14:08:46 2019 -0700

    Resolved pull conflict

[33mcommit 342e3d3ee94238305ebe8d444112ffbf40abb036[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri Jul 12 14:03:42 2019 -0700

    Snapshot on <sleec> before pull

[33mcommit 12d8dcdade54abba82a1f8b554b70ec3b92820fc[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Fri Jul 12 13:51:22 2019 -0700

    <sampling_effects.ipynb>, <sleec_input.py>:  Added <field_stop>
    option.  This demonstrates that adding a circular field stop fixes
    the artifacts with 4-fold symmetry seen in the contrast images.

[33mcommit 0a7da3df61301710b613557afcfccbda7f75068a[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Fri Jul 12 10:47:01 2019 -0700

    <sampling_effects.ipynb>:  Added scenario doubling focal-plane
    sampling (q=4 -> q=8).  No discernable improvement.  Also tried
    q=2, which noticeably degraded the contrast image.

[33mcommit d79fded5153ed99d05f0479763920d3dc0aa3d38[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Thu Jul 11 18:16:48 2019 -0700

    <sleec_input.py>:  Added <oversmpl> argument to <input_args()>
    <sampling_effects.ipynb>:  Adding double focal-plane sampling
    (..._2p_2f) scenario.  Paused at the end of the 7th cell.

[33mcommit d6846e24439fa870d8950213d3c8970cbd808c10[m
Merge: 8c443e6 4a3c291
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Thu Jul 11 18:06:08 2019 -0700

    Merge branch 'master' of https://github.jpl.nasa.gov/pinchen/sleec

[33mcommit 8c443e6fd51e02d58d9041e587c072bb09749794[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Thu Jul 11 18:02:58 2019 -0700

    Changed filname "pupil_pixels.ipynb" -> "sampling_effects.ipynb"
    Changed variable names in <sampling_effects.ipynb>:
        "..._4x" -> "..._2p"
        "..._1x" -> "..._1"

[33mcommit b39ec820b99a784d2f31107686e0aea1aee6cf24[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Thu Jul 11 17:29:29 2019 -0700

    <pupil_pixels.ipynb> last modified 2019-05-24:  Notebook to investigate relationship
    (if any) between pupil sampling density and artifaces in
    the coronagraph camera image.

[33mcommit 4a3c291c212503b8567bbba160e6b51935b56640[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Thu Jul 11 17:29:29 2019 -0700

    <pupil_pixels.ipynb>:  Notebook to investigate relationship
    (if any) between pupil sampling density and artifaces in
    the coronagraph camera image.

[33mcommit 430989b57f5bed2f03344b596a653fcb0adc299a[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Mon May 20 18:02:46 2019 -0700

    Snapshot

[33mcommit 3738a246c22017d9347118e773e0c37eb8317306[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Tue May 7 15:36:34 2019 -0700

    <sleec_input>: Improved <_try_mystic()> a bit.

[33mcommit 2db9151248bdb9f4b33fb88030bbc1dca48ee450[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Mon May 6 22:27:45 2019 -0700

    <sleec_input.py>: <_try_mystic> function is working.

[33mcommit 96f923bed3e4a81e7c5e6d34a1f011a95a036707[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Mon May 6 18:07:37 2019 -0700

    <sleec_input.py>: Modified <cost_ftn> (& <_test_cost_ftn>)
    functions to take FPM and Lyot-stop diameters as arguments.

[33mcommit f3a181464efa7072053644747fd36570f942badc[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Tue Apr 30 15:07:36 2019 -0700

    <sleec_input.py>: <cost_ftn> function seems to be working.
    See <_test_cost_ftn> function.

[33mcommit 22d0e756d155fd65d9bed2140fa2567add1c4181[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Fri Apr 26 17:18:07 2019 -0700

    <sleec_input.py>: Added <strehl> and <mean_contrast> functions.

[33mcommit a52b32f6d225bb8658e36499fe48de3c2d74739a[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Fri Apr 26 11:22:05 2019 -0700

    <sleec_optics.py>:  <Coronagraph> class.  Changed propagation
    output from list to dictionary.

[33mcommit 76fa8497a01cb4bd4fc4d0a242a01e170cd8c599[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Thu Apr 25 17:41:39 2019 -0700

    <sleec_input.py>: <ini_corngrph> function now takes an argument.
    Added <default_args> function.

[33mcommit 216014aa566c7653551bdc2eb02814b82aa5529b[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Thu Apr 25 14:24:08 2019 -0700

    <sleec_input.py>: Coded <lyot_thrghpt> function, passed
    <_test_lyot_thrghpt>.

[33mcommit ed078aae2b722db1f6e9aa8a18dd33e8fdd6dbc0[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Tue Apr 23 18:15:01 2019 -0700

    <sleec_input.py>: <fpm_thrghpt> function is working.  See
    <_test_fpm_thrghpt>.

[33mcommit 19b5883bf7bc8456335d62d4732b078017920882[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Mon Apr 22 20:55:48 2019 -0700

    <sleec_optics.py>:  Deleted <OptSystPlanesOut> class and absorbed its
    methods into <Coronagraph> class.  <Coronagraph> class:  added
    <load_elements> option to propagation methods.
    <sleec_input.py>:  Added <demo2> function to test the above changes.

[33mcommit d679974129d2679e5d4313a14ca7ad80bcc61d85[m
Merge: c36b900 810f01c
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Mon Apr 22 15:35:49 2019 -0700

    Merge branch 'dm'

[33mcommit 810f01c72cba3c2fdf3013247c81ceeb31aa5be9[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Fri Apr 19 14:26:00 2019 -0700

    <sleec_input.py>:  Revised <pupil_lens> function into <sphrcl_vmap>.
    <sleec_optics.py>:  After applying new <vmap> and updating <dm1>, need to
    update <optical_elements[1]>, same as before with updating <fpm>.

[33mcommit c36b900c1e7b9280649db17bf2199e92ec391b97[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Fri Apr 19 14:26:00 2019 -0700

    <sleec_input.py>:  Revised <pupil_lens> function into <sphrcl_vmap>.
    <sleec_optics.py>:  Need to debug why chaning vmap does not change DM
    surface.

[33mcommit bc33b58c8be4e8fc2a4be9c01319141334861d94[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Thu Apr 18 17:58:28 2019 -0700

    <sleec_input.py>.  Deleted <pupil_lens> function.  Coded
    <lens_wavefront>.  Need to test it.

[33mcommit 6d63adb00d24df37a8dd3bbaf19448d251cd8163[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Fri Apr 12 22:09:41 2019 -0700

    <sleec_input.py> <pupil_lens> function seems to work (qualitative
    check).

[33mcommit 2e6aee8c348f28966f047f20af205375b090e8ff[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Wed Apr 10 21:46:21 2019 -0700

    Propagation using physical units (instead of wvl=1, D=1, f=1) seems to
    working now, see <sleec_optics.py> <demo_coronagraph> function.

[33mcommit 86a363a29b8a8bec2a47db420457db43233f4978[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Tue Apr 9 22:27:41 2019 -0700

    Need to debug <sleec_optics.py> L571-572 (wavelength assignments
    in FraunhoferPropagator.

[33mcommit ee0f729468a2e33e5c2295e2178ef782e4d7dd36[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Tue Apr 9 17:39:11 2019 -0700

    Snapshot.  Need to test propagators w/ phys. units

[33mcommit 7578cb0081e2239e1680c3e5623e87be39c1fa4a[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Tue Apr 9 15:50:50 2019 -0700

    <sleec_input.py>:  Added <pupil_lens> function, untested.

[33mcommit c310328ab4501b9e922b707272e3c4fc4908191f[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Thu Apr 4 18:24:58 2019 -0700

    Renamed <hcipy> subdirectory to <hcipy_local>.
    <sleec_input.py>:  started coding throughput functions (incomplete).

[33mcommit 522c74d8102899bdb4cb92a80964b610eba72a6c[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Apr 2 19:05:34 2019 -0700

    <sleec_input.py> and <sleec_optics.py>:  Code clean-up.

[33mcommit d9d19be4014e37baf4412728eed5cc4175c40ae6[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Apr 2 19:01:34 2019 -0700

    <sleec_input.py>:  Fixed bug of <update_fpm> function not taking effect.  Had
    to load the <.fpm> attribute for the update to fully take effect (due to
    @property lazy loading).

[33mcommit a517a0ec667e2d646d304f708938b4f53d79ce39[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Mon Apr 1 17:46:17 2019 -0700

    <sleec_input.py>:  Added <update_fpm> function.  <sleec_optic.py>:  Need to
    modify <OptclSysPlanesOut> class to reset <wavefronts> attribute after
    executing <update_fpm>.

[33mcommit c0517d0ec4aec4c3b9f17ddab1281cbb253ba4b3[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Mar 27 17:50:47 2019 -0700

    Renamed <coro_sleec.py> -> <sleec_optics.py>

[33mcommit bca34512eae989afc826ffb149df34f288d700d2[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Mar 27 17:42:28 2019 -0700

    <sleec_input.py>:  Created <input_args> function containing line that
    populate the <pars> dictionary (with fields and DM voltage-maps at optical
    planes.
    Moved <coro_sleec.py><Coronagraph.load_args> method to become
    <sleec_input.py><ini_corngrph> function.

[33mcommit 53436e63a6768af8c6c9034d5863b52bcbe0381c[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Mar 27 12:04:35 2019 -0700

    <sleec_input.py>:  Started coding to enable manipulation/mutation of the
    coronagraph object.

[33mcommit 0b6ca38e3c781732f33b55d578be639a8bfc114d[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri Mar 22 18:14:55 2019 -0700

    <coro_sleec.py>, <Coronagraph> class:  Replaced <make_elements> method with
    <@property> attributes so these parameters can be updated on the fly.

[33mcommit d4e9a3e1c4827d4747e1115d8ccf84401fe85558[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Mar 21 17:49:36 2019 -0700

    <coro_sleec.py>: Replaced <_test()> function with <demo_coronagraph()>.

[33mcommit 4cbbd4a25cd2a5380b38c8b584e07f2ed80837cb[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Mar 21 17:00:48 2019 -0700

    <coro_sleec.py>, <Coronagraph> class: Removed extraneous attributes.

[33mcommit e4176f94a9b99b99a2ac801357f51dab2c08ebb6[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Mar 21 16:54:21 2019 -0700

    <coro_sleec.py>: New class <Coronagraph> now working.  <sleec_input.py> OK
    too.

[33mcommit 02c5a517737e914a15006ce3638db85e0bcd4fe3[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Mar 20 18:02:11 2019 -0700

    <coro_sleec.py>:  Snapshot.

[33mcommit 3fde3feb4e746b19247010d577e092df40cb410f[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Mar 20 15:53:13 2019 -0700

    <coro_sleec.py>:  <TransmMask> class, <input_grid> and <wavelength> arguments
    are not necessary.  Removed them.

[33mcommit 326072edebe9b5349daaeb2bd2722f7bff6b783d[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Mar 20 15:44:14 2019 -0700

    <coro_sleec.py>:  Snapshot before trying to get rid of some extraneous(?)
    variables.

[33mcommit d222f0dad5f4322ee027332258ecb0841dddeb3f[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Mar 19 22:07:59 2019 -0700

    Snapshot

[33mcommit 0b93a6b42ba74029ec21b9fb7db7dea6e4bbfc11[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Mon Mar 18 22:13:28 2019 -0700

    Started coding <sleec_input.py>.  Need to complete it, and then remove
    redundant code in <coro_sleec.py><Coronagraph> class.

[33mcommit f47f5002bd243e941806b91b4c4731d2ad1575ee[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri Mar 15 18:59:48 2019 -0700

    Just migrated to newly replenished computer.
    <coro_sleec.py>:  Started coding <Coronagraph> class based on <_test()>
    function.

[33mcommit 39903a7a5ea119eda25a9dce73168dafda818290[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Mar 7 15:00:43 2019 -0800

    <coro_sleec.py>:  <_test> function now returns list of names of optical planes

[33mcommit 311276ca19f41bdf1de3c9d0f42d5de44cce36b9[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Mar 6 17:57:18 2019 -0800

    <coro_sleec.py>:  Added <OptclSysPlanesOut> class.  It is the same as
    hcipy.OpticalSystem, but it outputs the propagated wavefront from every
    optical element.  <_test()> method now uses this new class to output all
    wavefronts.

[33mcommit 57ac32e78e979fee09e265136c67996267ecf60e[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri Feb 8 17:14:51 2019 -0800

    <coro_sleec.py>:  Tidied up code a bit.  <_test()> now returns variables
    instead of plotting results.

[33mcommit 82c7431d33092c4896b48567fed81c19dc473de7[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Feb 7 12:52:32 2019 -0800

    Started <verify_coro_sleec.ipynb> to verify accuracy of propagation
    computations.

[33mcommit 6e5a54582c3cc812677b7ee9bb8b9ed8e54562f5[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Feb 6 18:13:54 2019 -0800

    Inserted emporary code block (starting L467) to examine influence-function
    overlap effect.

[33mcommit 4482c4bbd93f0437662334ade2d924550d365fb6[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Mon Feb 4 18:35:05 2019 -0800

    <coro_sleec.py>:  Fixed some tab vs. space formatting issues.

[33mcommit c44d12ae75d397e39528bdd753da17358f86a708[m
Merge: c38a73b c54804d
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri Feb 1 17:55:10 2019 -0800

    Merge branch 'master' of https://github.jpl.nasa.gov/pinchen/sleec

[33mcommit c38a73bb9c6a80955a7cfc16e5e0fc6e3ca01fd4[m
Merge: 659fec8 16cb80a
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri Feb 1 17:53:30 2019 -0800

    Merge branch 'master' of https://github.jpl.nasa.gov/pinchen/sleec

[33mcommit c54804d0f845a83e1b490bf5c16f47fb0ae04e29[m
Merge: 659fec8 16cb80a
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri Feb 1 17:53:30 2019 -0800

    Merge branch 'master' of https://github.jpl.nasa.gov/pinchen/sleec

[33mcommit 659fec88e1d0e3a2e42135444f7e02912233e5cf[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri Feb 1 17:51:04 2019 -0800

    <coro_sleec.py> _test method:  Modified to run optical propagation using
    sparse-matrix influence functions.  Trid flat DMs and random actuation on DM1.
    At least qualitatively, propagation appeared to work fine.

[33mcommit eed4e3be9d597c2af950821a94de7dd9fec1d79b[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri Feb 1 16:34:56 2019 -0800

    <coro_sleec> InflncFtn class:  Revised <generate> (renamed to <genrt_modes>)
    and <load_results> (renamed to <load_modes>) methods to
    keep their outputs as sparse matrices.  Added <make_surfc> method to construct
    DM surface using sparse-matrix influence functions.

[33mcommit 16cb80a19c700d3a7e7abca43511f562ecf205ff[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri Feb 1 16:34:56 2019 -0800

    <coro_sleec> InflncFtn class:  Revised <generate>, <load_results> methods to
    keep their outputs as sparse matrices.  Added <make_surfc> method to construct
    DM surface using sparse-matrix influence functions.

[33mcommit 027f08793fab800b40dc356167c2bfa8100d4513[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Jan 31 18:10:57 2019 -0800

    <coro_sleec.py> InflncFtn.load_results():  convert sparse matrix to np.array.

[33mcommit eb90ec3bd66f2143d2728479b0009cfacd80fc7c[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Jan 31 17:58:01 2019 -0800

    <coro_sleec.py>  InflncFtn class:  <generate> method now saves a dictionary
    of the main output (sparse CSR matrix of influence functions) and relevant
    parameters.  Coded <load_results> method.  It loads the pickled resultes
    from <generate> and constructs a hcipy.ModeBasis of influence functions.  The
    fields of the ModeBasis are sparse matrices.  I expect this to fail at some
    point in optical progagations.  Using the <.toarray()> method will fix it.

[33mcommit d3389eea86a4fba600fc6e96d3a921ef971d3b0c[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Jan 30 15:30:39 2019 -0800

    <coro_sleec.py> InflncFtn class:  Implemented <make_intrpltr> method, which
    emulates <hcipy.make_linear_interpolator_separated> method with added option
    of linear or nearest interpolation.
    Test with DM48.4, 10x10 super-influence function (superinf-20180730-10x10.fits)
    as input, 481-pixel pupil grid:
    - The nearest option exactly reproduces the original (input) super-influence
    function.
    - The linear method reproduces to within 3.8E-15.
    
    Load time of the saved NPZ file (43 MB) = 1.7 s

[33mcommit 1cdb41b5e4681cc456606788744a873ad844049e[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Jan 30 11:19:34 2019 -0800

    <coro_sleec.py> InflncFtn class:  Made adjustments to how grids are defines
    so that coordinates in <self.acttr.grid> and <self.pupil_grid> overlap.

[33mcommit 52bf8a01c7e32948a4b0460c466607b309fc3afc[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Mon Jan 28 17:55:27 2019 -0800

    <coro_sleec.py> InflncFtn class:  Minor adjustments to decorators, attributes,
    and arguments.

[33mcommit 1c21921a70816873f083f7ae63ee46a38721a4c6[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Mon Jan 28 10:51:54 2019 -0800

    <coro_sleec.py> InflncFtn.acttr_grid:  Interpolation scaling now looks
    OK, with single-actuator influence function spanning about 6 pxls HWHM.

[33mcommit 658ef0a2f02060084e9577754a309385a5257130[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri Jan 25 22:14:35 2019 -0800

    <coro_sleec.py> InflncFtn.acttr_grid: Updated 'Notes' docstring.

[33mcommit 6f7c84b41f2bae2f4fd9e0744d676ea4f0f35764[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri Jan 25 18:04:21 2019 -0800

    <coro_sleec.py> InflncFtn class:  switched to using 10x10 super-influence
    function.

[33mcommit c237f7bb21de677cd223601517589d893b96f27a[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Jan 24 18:23:59 2019 -0800

    <coro_sleec.py> InflncFtn class:  Saving using sparse matrix works now.

[33mcommit fd894ccb235908d1567f3c12d1ead968ab601075[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Jan 24 14:24:10 2019 -0800

    <coro_sleec.py> InflncFtn class:  Tweaked definitions of <scale_row> and
    <scale_column> arguments.  Moved <pupil_grid> argument to class interface.

[33mcommit d882cb2a5ca74cb79dea3e50246e6708c12b4c34[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Jan 23 16:25:27 2019 -0800

    <coro_sleec.py> InflncFtn.acctr_grid:  Figured out how to properly set values
    for  <scale_row>, <scale_column> arguments.

[33mcommit eff91adb265beed13e40164384510dd51ff930b7[m
Merge: e397de5 6733bf0
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Jan 22 16:51:47 2019 -0800

    Merge branch 'master' of https://github.jpl.nasa.gov/pinchen/sleec

[33mcommit e397de5e8b1bd40b4b62667b87c4737f9aa2e615[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Jan 22 16:47:31 2019 -0800

    <coro_sleec.py>:  After debugging, <InflncFtn> class runs w/o runtime error
    now.  Still need to verify the class and implement sparse matrices.

[33mcommit 6733bf0c3e7753e5bfea152b716880b8393e987f[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Jan 22 16:47:31 2019 -0800

    <coro_sleec.py>:  After some debugging, no runtime error now.  Still need to
    verify and implement sparse matrices.

[33mcommit 88c85b1199377256d045878fbd065bbd287462b0[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri Jan 18 17:48:19 2019 -0800

    <coro_sleec.py>
    Coded <InflncFtn> class.  Need to test run.

[33mcommit 5bad3f1a60fcd0588ab44c52576e8523a74897b0[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Jan 17 18:24:52 2019 -0800

    <coro_sleec.py>:  Working on using sparse matrices for storing DM influence
    functions to reduce data-loading time.  Need to finish coding <InflncFtn>
    class.

[33mcommit be73cd39ba9708a6cc1984597c05f0b751197ff9[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Jan 16 15:09:28 2019 -0800

    <coro_sleec.py>:  Removed <SleecCoronagraph> class and <_test_class> function.

[33mcommit 67154d82e6e6441293dc9a1ddc9b9e805f7f5a55[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Jan 16 12:21:55 2019 -0800

    Snapshot

[33mcommit 5781fa54a59464120cb86500a5757101c40d8846[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Dec 20 16:35:35 2018 -0800

    <learn_mystic.py>
        Make small changes to improve plotting.

[33mcommit 6a753f29f0dc787a4c0ae7dd31a0524bfe88ba32[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Dec 20 14:44:22 2018 -0800

    <learn_mystic.py>
    Made improvements.

[33mcommit 7870c77cffa837923b34b5c325fe5b7ab1ade2c3[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Dec 20 11:00:55 2018 -0800

    <learn_mystic.py>
    Revised code to use <CostFactory2> of mystic <AbstractModel>

[33mcommit 7739d7aa6b868dca8fd0ad64988d3f510940cf8f[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Dec 19 18:28:56 2018 -0800

    <learn_mystic.py> appears to be working.

[33mcommit 434a06557ebbb12ad4c9c1fee27727ccb3b09d33[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Dec 18 18:13:30 2018 -0800

    <learn_mystic.py>
    Wrote trial mystic model.  Not tested yet.

[33mcommit 92cfd0ebaa7ea0dba772d9511fde8d2db836f633[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri Dec 14 17:52:59 2018 -0800

    <learn_mystic.py>:  Incomplete

[33mcommit 4b3515d1ad58bdb89860adf3391cb5482e058135[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Dec 11 15:01:38 2018 -0800

    <mystic_example06.py>:  Removed extraneous lines.

[33mcommit d266f995f60bd4f32bde1ca6b0fd48f63179700a[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Dec 11 14:29:28 2018 -0800

    Copied <example06.py> from <mystic>.  Replaced <pylab> in the script by
    <matplotlib.pyplot>. Runs w/o error.

[33mcommit a84696df2458882050148e9333e81f59c76a0963[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri Dec 7 17:25:05 2018 -0800

    Rebuilt Sphinx docs

[33mcommit 937e65bd0991ab68f69eb3ba00dbc68253eddd0a[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri Dec 7 16:47:47 2018 -0800

    After learning about serious pitfalls of git modules, decided to remove the
    <hcipy> module and copy it back as an integral part of <sleec> code.  Did the
    following:
    1. Created an archive copy of <sleec> named <sleec_backup_20181207>.
    2. In <sleec/source/> directory, issued command:
        $ git rm hcipy
    3. Removed <modules> subdirectory under <.git> directory and .gitmodules (in
       <sleec> directory.
    4. Copied back <hcipy>.
    5. Removed all .git and .gitignore files in <hcipy>.  These were the only
       .git-type files in <hcipy>.
    6. $ git add -A :/

[33mcommit be23c1f5a3b7ee611599d833bd055a400aad331f[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Dec 6 11:43:45 2018 -0800

    <coro_sleec.py>
    Current version generated JPEG images to show effects of DM1 and DM2 by
    putting sinusoidal on one vs. the other.

[33mcommit 969af81e58d07a8aa3722838b898bef466192778[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Dec 6 09:52:09 2018 -0800

    <coro_sleec.py>:  <_test_instance>  pickling DM influence-functions mode
    basis.  Reduces run time 15x vs. generating on the fly.

[33mcommit 4add0abe280a8e0816c4cb956895b6924df20069[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Dec 5 18:30:38 2018 -0800

    <coro_sleec.py>
        Added two DMs to <_test_instance()> method.  Seems to be working:  put
        random settings on DM1 actuators, speckles show up in final image.

[33mcommit 6cf89905ab4fc71c2a7c127a6b17da2671384d61[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Dec 5 16:49:23 2018 -0800

    <coro_sleec.py>
        <_test_instance> is working --> Can model coronagraph as an
        <hcipy.OpticalSystem>.

[33mcommit 75c5e8661c517312272abfaccaf7939a5a956203[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Dec 5 09:50:45 2018 -0800

    Disabled sphinx autosummary.  Fixed docstring bugs in <coro_sleec.py>.

[33mcommit 594e77fd61bc73cfc1069fdbdbb4ad4182db8b28[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Dec 5 09:38:00 2018 -0800

    Rebuilt sphinx docs to try autosummary extension.  Doesn't seem to make a
    difference.

[33mcommit 1f7565cad9774196ce764e40f91219f52155b3e7[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Dec 4 23:29:53 2018 -0800

    Rebuilt sphinx documentation.
    <coro_sleec.py>:  Updated <_test()>.  Determined dark outer ring in the
    intensity field before the Lyot plane is an artifact of small field of view
    for the focal-plane grid.

[33mcommit 858e23d4610d8b7b3ece73ab67f729a755bcaa05[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Dec 4 16:37:37 2018 -0800

    <coro_sleec.py>:  Renamed <FocalPlaneMask> to <TransmMask>.  Edited it to
    emulate <hci.ApodizerMonochromatic>, with added <input_grid> attribute.

[33mcommit 3fe847a8f82e65efd66d1b36f94959974118d998[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Dec 4 14:39:27 2018 -0800

    Back story:  Installed <sphinx> to generate documentation pages.  Installed
    <sphinx> in <sleec/doc> subdirectory.  Created <sleec/source> subdirectory as
    code/package directory.
    Changes since last commit:  Had to mode <hcipy> package into <source> as well
    to make <sphinx/autooc> work.  Otherwise, <import hcipy> statements in my
    code cause errors.  Had to edit path specifications in the following files to
    make git work again:
    ~/sleec/.gitmodules
    ~/sleec/source/hcipy/.git
    ~/sleec/.git/modules/hcipy/config

[33mcommit d356ee8d0a5af321998593cdf954ef706fbc1f56[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Tue Dec 4 13:44:38 2018 -0800

    Snapshot before messing further with directory structure to debug sphinx.

[33mcommit 88911d63186c8f05191621194321a8ce7e0361c0[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri Nov 30 22:22:04 2018 -0800

    <hci.make_xinetics_influence_functions()> now works in <2dm_coronagraph.ipynb>

[33mcommit 4c5a45a6699254e520fabc3f6e52e0e02994ab00[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri Nov 30 21:41:46 2018 -0800

    Snapshot

[33mcommit ecb3bbececd96bfdc0e3da870bd883dcf5e79ba2[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri Nov 30 17:52:00 2018 -0800

    Added <hcipy> submodule

[33mcommit be05c68a9dd577fa53aae15d3fda5c1806266abc[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Nov 28 17:06:57 2018 -0800

    Same as before, but saved ipynb checkpoint.

[33mcommit e9cba44a598b58f058f87bd0b6aa13d8e4196c53[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Wed Nov 28 17:03:39 2018 -0800

    <2dm_coronagraph.ipynb>:  Working without runtime error.

[33mcommit 0133c03d7cfa5bdbfae0ad925ee462770a54483f[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Thu Nov 22 14:03:09 2018 -0800

    Edited README.md following the style in <README-Template.md>, downloaded from
    <https://gist.github.com/PurpleBooth/109311bb0361f32d87a2>.

[33mcommit c4f76459aa1c6606d15ac62ce282e448203e5cc9[m
Author: Pin Chen <pin.chen@jpl.nasa.gov>
Date:   Fri Nov 16 15:04:08 2018 -0800

    Created <learn_hcipy.ipynb> as my workbook in learning HCIPy.

[33mcommit f8b3e4df26210d0344687bc6ab79ca5b9dc364b0[m
Author: Pin Chen <Pin.Chen@jpl.nasa.gov>
Date:   Fri Nov 16 12:39:18 2018 -0800

    Initial commit
