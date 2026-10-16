def read_multiline_fasta(input_fasta: str) -> list[tuple[str, str]]:
    """
    Читает многострочный FASTA-файл и объединяет разбитые последовательности.

    Arguments:
        input_fasta (str): Путь к входному FASTA-файлу.

    Returns:
        list[tuple[str, str]]: Список кортежей, где первый элемент — заголовок (>ID),
                               второй элемент — полная однострочная последовательность.
    """
    results: list[tuple[str, str]] = []
    current_id: str = ""
    current_seq_chunks: list[str] = []

    with open(input_fasta, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue  # Пропускаем пустые строки, если они есть

            # Если строка начинается с '>', значит это новый заголовок
            if line.startswith(">"):
                # Если у нас уже накопилась предыдущая последовательность, сохраняем её
                if current_id:
                    full_seq = "".join(current_seq_chunks)
                    results.append((current_id, full_seq))

                # Сбрасываем буфер для нового рида
                current_id = line
                current_seq_chunks = []
            else:
                # Если это строка с нуклеотидами/аминокислотами, добавляем её в кусочки
                current_seq_chunks.append(line)

        # Сохраняем самый последний рид после выхода из цикла
        if current_id:
            full_seq = "".join(current_seq_chunks)
            results.append((current_id, full_seq))

    return results
