namespace {
void fillMissingValuesNearest(const atlas::FieldSet & sourceFieldSet,
                              atlas::FieldSet & targetFieldSet,
                              const oops::Variables & vars,
                              const atlas::FunctionSpace & sourceFs,
                              const atlas::FunctionSpace & targetFs) {
  if (vars.size() == 0) {
    oops::Log::trace() << classname() << "fillMissingValuesNearest: no variables to process" << std::endl;
    return;
  }

  const auto src_lonlat = atlas::array::make_view<double, 2>(sourceFs.lonlat());
  const auto src_ghost = atlas::array::make_view<int, 1>(sourceFs.ghost());
  std::vector<double> lons;
  std::vector<double> lats;
  std::vector<atlas::idx_t> indices;
  lons.reserve(src_lonlat.shape(0));
  lats.reserve(src_lonlat.shape(0));
  indices.reserve(src_lonlat.shape(0));
  for (atlas::idx_t jj = 0; jj < src_lonlat.shape(0); ++jj) {
    if (src_ghost(jj) == 0) {
      lons.push_back(src_lonlat(jj, 0));
      lats.push_back(src_lonlat(jj, 1));
      indices.push_back(jj);
    }
  }
  if (indices.empty()) {
    oops::Log::trace() << classname() << "fillMissingValuesNearest: no owned source points" << std::endl;
    return;
  }

  const atlas::Geometry earth(atlas::util::Earth::radius());
  atlas::util::IndexKDTree2D tree(earth);
  tree.build(lons, lats, indices);

  const auto tgt_lonlat = atlas::array::make_view<double, 2>(targetFs.lonlat());
  const auto tgt_ghost = atlas::array::make_view<int, 1>(targetFs.ghost());
  const double missing = util::missingValue<double>();


  for (const auto & var : vars) {
    if (!targetFieldSet.has(var.name()) || !sourceFieldSet.has(var.name())) {
      oops::Log::trace() << classname() << "fillMissingValuesNearest: skipping var (missing in fset)" <<var.name()<< std::endl;
      continue;
    }
    auto tgt_view = atlas::array::make_view<double, 2>(targetFieldSet[var.name()]);
    const auto src_view = atlas::array::make_view<double, 2>(
        sourceFieldSet.field(var.name()));

    std::size_t missing_before = 0;
    std::size_t missing_after = 0;
    std::size_t filled = 0;
    std::size_t logged = 0;
    const bool log_values = (var.name() == "air_pressure_thickness");
    const std::size_t log_limit = 20;
    const double small_value_threshold = 1.0e-6;
    std::size_t small_after_fill = 0;
    std::size_t small_logged = 0;

    for (atlas::idx_t jloc = 0; jloc < tgt_view.shape(0); ++jloc) {
      if (tgt_ghost(jloc) != 0) {
        continue;
      }
      bool has_missing = false;
      for (atlas::idx_t jlev = 0; jlev < tgt_view.shape(1); ++jlev) {
        if (tgt_view(jloc, jlev) == missing) {
          has_missing = true;
          ++missing_before;
          break;
        }
      }
      if (!has_missing) {
        continue;
      }

      
    }

    for (atlas::idx_t jloc = 0; jloc < tgt_view.shape(0); ++jloc) {
      if (tgt_ghost(jloc) != 0) {
        continue;
      }
    }

  }
