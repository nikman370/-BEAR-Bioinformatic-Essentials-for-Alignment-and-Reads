def filter_fastq(
    input_path: str,
    output_path: str,
    gc_bounds: tuple[float, float] | float = (0, 100),
    length_bounds: tuple[int, int] | int = (0, 2**32),
    quality_threshold: float = 0,
) -> None:
    """
    Главная функция для фильтрации FASTQ-файлов.
    Читает input_path, проверяет каждый рид и сохраняет прошедшие в output_path.
    """
