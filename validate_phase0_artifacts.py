from __future__ import annotations

import ast
import json
from pathlib import Path

ROOT = Path(__file__).parent
expected = [
    '00_setup_aoi_and_config.ipynb',
    '01a_gee_grace.ipynb', '01b_gee_smap.ipynb',
    *[f'01{letter}_gee_{name}.ipynb' for letter, name in [
        ('d', 'tropomi_sif'), ('e', 'era5_land'), ('f', 'gridmet'),
        ('g', 'modis_ndvi_evi'), ('h', 'modis_landcover'), ('i', 'srtm_dem'),
        ('j', 'usda_cdl'), ('k', 'gfsad'), ('l', 'aridity_index'),
        ('m', 'worldcover'), ('n', 'landsat'), ('o', 'sentinel2'),
        ('p', 'pml_v2'), ('q', 'ssebop'), ('r', 'soilgrids'), ('s', 'nasadem')]],
    *[f'02{letter}_indep_{name}.ipynb' for letter, name in [
        ('a', 'ma_wtd_raster'), ('b', 'usgs_wells'), ('c', 'state_wells'),
        ('d', 'ameriflux'), ('e', 'sapfluxnet'), ('f', 'gleam4'),
        ('g', 'gosif'), ('h', 'lanid'), ('i', 'nass_stats'), ('j', 'ecostress')]],
    '98_monitor_gdrive_exports.ipynb', '99_run_all_and_validate.ipynb',
]
errors: list[str] = []
for name in expected:
    path = ROOT / name
    if not path.exists():
        errors.append(f'missing: {name}')
        continue
    try:
        document = json.loads(path.read_text(encoding='utf-8'))
    except Exception as exc:
        errors.append(f'{name}: json: {exc}')
        continue
    cells = document.get('cells', [])
    if len(cells) != 12:
        errors.append(f'{name}: expected 12 cells, got {len(cells)}')
    if not cells or cells[0].get('cell_type') != 'markdown':
        errors.append(f'{name}: first cell is not markdown')
    for index, cell in enumerate(cells, start=1):
        if 'language' not in cell.get('metadata', {}):
            errors.append(f'{name}: cell {index} lacks metadata.language')
        if cell.get('cell_type') == 'code':
            source = ''.join(cell.get('source', []))
            source = '\n'.join(line for line in source.splitlines() if not line.startswith(('%', '!')))
            try:
                ast.parse(source, filename=f'{name}:cell{index}')
            except SyntaxError as exc:
                errors.append(f'{name}: cell {index}: {exc}')
    if name.startswith('01'):
        text = path.read_text(encoding='utf-8')
        for required in ['ee.batch.Export.image.toDrive', 'BatchManager', 'MASTER_CRS', 'MASTER_SCALE', 'MAX_CONCURRENT = 20']:
            if required not in text:
                errors.append(f'{name}: missing {required}')
        for forbidden in ['google.colab', 'drive.mount', '/content/drive', 'write_exported_geotiff_as_cf']:
            if forbidden in text:
                errors.append(f'{name}: forbidden local-mode token {forbidden}')

config = json.loads((ROOT / 'config.yaml').read_text(encoding='utf-8'))
for dataset in config.get('datasets', {}).values():
    if not dataset.get('asset_id'):
        errors.append('config: dataset missing asset_id')
for dataset in config.get('independent_datasets', {}).values():
    if not dataset.get('name'):
        errors.append('config: independent dataset missing name')

print(f'expected_notebooks={len(expected)}')
print(f'validated_notebooks={len([name for name in expected if (ROOT / name).exists()])}')
print(f'validation_errors={len(errors)}')
for error in errors[:20]:
    print(error)
raise SystemExit(1 if errors else 0)
